# scripts/test_wrappers/filter_species_set_6sp_test.R ---------------------------
#
# Stage 03: Policy filtering for a small species set (pre-rasterisation)
#
# Purpose:
#   Run Stage 03 over a short species list using a single, explicitly-defined policy.
#   Stage 03 is expected to be re-run frequently as thresholds change; outputs are
#   overwritten by default to avoid stale "filtered" datasets persisting.
#
# Inputs:
#   data/processed/02_qc_flagged/<slug>/occ_<slug>__qc_flagged.(parquet|rds)
#
# Outputs:
#   data/processed/03_filtered/<slug>/occ_<slug>__filtered.(parquet|rds)
#   data/processed/03_filtered/_runlog_03_filtered.csv              (optional)
#
# Policy update (Feb 2026):
#   Include specimen/museum-type records provided they have usable coordinates
#   and are within the normal date window (post-2000 for now).
#   For GBIF this is implemented via basisOfRecord inclusion + a simple provenance gate
#   for specimen/material-sample rows.
#
# ------------------------------------------------------------------------------

# ---- Find repo root (works from any scripts/ subfolder) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file) || !nzchar(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run via source('.../scripts/.../filter_species_set_6sp_test.R') from a file, not copy/paste."
  )
}

script_dir <- dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE))

find_repo_root <- function(start_dir) {
  markers <- c(".git", "R", "data", "InfluentialSpecies.Rproj", "DESCRIPTION")
  d <- start_dir
  for (i in 1:15) {
    if (any(file.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root <- find_repo_root(script_dir)

# ---- Load the Stage 03 engine -------------------------------------------------
engine_fn <- file.path(repo_root, "R", "filter_occurrences.R")
if (!file.exists(engine_fn)) {
  stop("Can't find Stage 03 engine at: ", engine_fn)
}
source(engine_fn)

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# ---- Species list -------------------------------------------------------------
species_names <- c(
  "Myrmica sabuleti",
  "Myrmica scabrinodis",
  "Andrena fulva",
  "Sorex araneus",
  "Leptothorax acervorum",
  "Emberiza schoeniclus"
)

# ---- Stage 03 policy ----------------------------------------------------------
# Stage 03 is policy, so every numerical decision lives here.
# The policy_id is written to the runlog and should change if/when thresholds change.
#
# Notes on basisOfRecord:
#   NBN often has basisOfRecord == NA, so we do NOT apply a global allowed_basis_of_record
#   gate here (that would accidentally drop large amounts of NBN).
#   Instead we enforce GBIF-specific basis/provenance rules via extra_drop_rules.

policy <- list(
  # Unique label for this exact set of filtering rules.
  # Change this whenever you change ANY threshold/switch (e.g. uncertainty 1 km -> 5 km),
  # so run logs and outputs can be traced back to the policy that produced them.
  #
  # This policy:
  #   - keeps post-2000 records (for now) with <=1 km uncertainty (where known)
  #   - includes in-situ observations AND specimen/material-sample records (GBIF)
  #   - requires specimen/material-sample records to have real provenance fields populated
  policy_id = "baseline_2000_unc1km_obs_plus_specimen_prov_gbif",
  
  # Keep only selected sources (NULL keeps all).
  keep_sources = NULL,              # e.g. c("GBIF", "NBN")
  
  # ---- Structural drops (usually TRUE) ----------------------------------------
  drop_missing_coords = TRUE,
  drop_coords_out_of_range = TRUE,
  drop_future_date = TRUE,
  
  # ---- Date completeness ------------------------------------------------------
  drop_missing_date = FALSE,
  
  # ---- Date window (applies to event_day; may fallback to year if allowed) ----
  min_date = "2000-01-01",          # set NA to disable
  max_date = NA,                    # set NA to disable
  allow_year_only = TRUE,
  require_event_day = FALSE,
  min_year = NA,
  max_year = NA,
  
  # ---- Coordinate uncertainty -------------------------------------------------
  max_coord_uncertainty_m = 1000,   # 1 km; set 5000 for 5 km, 10000 for 10 km, etc.
  uncertainty_missing_action = "keep",
  
  # ---- GBIF issues handling ---------------------------------------------------
  gbif_issues_mode = "ignore",
  issues_blacklist = c(),
  
  # ---- Licence handling -------------------------------------------------------
  drop_unexpected_licence = FALSE,
  
  # ---- basisOfRecord handling (global; see note above) ------------------------
  allowed_basis_of_record = NULL,
  drop_basis_of_record = NULL,
  
  # ---- Taxon rank gating (optional; applies if taxonRank exists) --------------
  allowed_taxon_rank = NULL,
  
  # ---- NBN certainty gating (optional; applies if source=="NBN" and column exists)
  nbn_certainty_col = NULL,
  nbn_allowed_certainty = NULL,
  
  # ---- Extra rules (advanced) -------------------------------------------------
  # Named list of functions(dt) -> logical drop vector (TRUE means drop).
  extra_drop_rules = list(
    
    # Rule 1: GBIF basisOfRecord inclusion gate.
    #
    # Keep only:
    #   - HUMAN_OBSERVATION
    #   - OBSERVATION
    #   - MACHINE_OBSERVATION
    #   - PRESERVED_SPECIMEN
    #   - MATERIAL_SAMPLE
    #
    # Drop everything else for GBIF (including missing basisOfRecord for GBIF).
    drop_gbif_basis_not_in_allowed_set = function(dt) {
      if (!("source" %in% names(dt)) || !("basisOfRecord" %in% names(dt))) {
        return(rep(FALSE, nrow(dt)))
      }
      
      is_gbif <- dt$source == "GBIF"
      bor <- toupper(trimws(as.character(dt$basisOfRecord)))
      
      allowed <- c(
        "HUMAN_OBSERVATION",
        "OBSERVATION",
        "MACHINE_OBSERVATION",
        "PRESERVED_SPECIMEN",
        "MATERIAL_SAMPLE"
        # If we later decide these should count too, add here:
        # "LIVING_SPECIMEN"
      )
      
      is_gbif & !(bor %in% allowed)
    },
    
    # Rule 2: Specimen/material-sample provenance requirement (GBIF only).
    #
    # We only enforce this for:
    #   - PRESERVED_SPECIMEN
    #   - MATERIAL_SAMPLE
    #
    # A record passes if it has at least one "real provenance" identifier present.
    # This is intentionally light-touch; on our merged schema Andrena fulva had 100%
    # provenance coverage across these fields upstream.
    drop_gbif_specimen_missing_provenance = function(dt) {
      if (!("source" %in% names(dt)) || !("basisOfRecord" %in% names(dt))) {
        return(rep(FALSE, nrow(dt)))
      }
      
      is_gbif <- dt$source == "GBIF"
      
      bor <- toupper(trimws(as.character(dt$basisOfRecord)))
      is_spec <- bor %in% c("PRESERVED_SPECIMEN", "MATERIAL_SAMPLE")
      
      has_val <- function(col) {
        if (!(col %in% names(dt))) return(rep(FALSE, nrow(dt)))
        x <- trimws(as.character(dt[[col]]))
        !is.na(x) & nzchar(x)
      }
      
      # Use the provenance fields that actually exist in our schema.
      prov_ok <- (
        has_val("occurrenceID") |
          has_val("datasetKey") |
          has_val("institutionCode") |
          has_val("collectionCode") |
          has_val("datasetName")
      )
      
      is_gbif & is_spec & !prov_ok
    }
  )
)

# ==============================================================================
# RUN SETTINGS
# ==============================================================================

in_root  <- file.path(repo_root, "data", "processed", "02_qc_flagged")
out_root <- file.path(repo_root, "data", "processed", "03_filtered")

overwrite     <- TRUE
write_parquet <- TRUE
write_rds     <- FALSE
write_runlog <- TRUE

continue_on_error <- TRUE
verbose <- TRUE

# ==============================================================================
# RUN
# ==============================================================================

stage03_filter_occurrences(
  species_names = species_names,
  policy = policy,
  in_root = in_root,
  out_root = out_root,
  overwrite = overwrite,
  write_parquet = write_parquet,
  write_rds = write_rds,
  write_runlog = write_runlog,
  continue_on_error = continue_on_error,
  verbose = verbose
)

# ==============================================================================
# Quick Stage 03 “before vs after” + drop-reason summary (prints to console)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

runlog_path <- file.path(out_root, "_runlog_03_filtered.csv")

if (!file.exists(runlog_path)) {
  message("[Stage 03 summary] No runlog found at: ", runlog_path)
} else {
  
  hdr <- strsplit(readLines(runlog_path, n = 1, warn = FALSE), ",", fixed = TRUE)[[1]]
  
  lg <- fread(
    runlog_path,
    sep = ",",
    header = TRUE,
    fill = TRUE,
    quote = "\"",
    na.strings = c("", "NA")
  )
  
  if (ncol(lg) > length(hdr)) {
    setnames(lg, c(hdr, paste0("extra_col_", seq_len(ncol(lg) - length(hdr)))))
  }
  
  if (!("timestamp_utc" %in% names(lg))) {
    message(
      "[Stage 03 summary] Runlog parsed unexpectedly (timestamp_utc missing).\n",
      "If this persists, delete and regenerate the runlog:\n",
      "  file.remove(runlog_path)\n",
      "Then re-run Stage 03."
    )
  } else {
    
    lg[, timestamp_utc_parsed := as.POSIXct(
      timestamp_utc,
      format = "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    )]
    
    slugify_local <- function(x) {
      s <- gsub("[^a-z0-9]+", "_", tolower(x))
      gsub("^_+|_+$", "", s)
    }
    slugs <- vapply(species_names, slugify_local, character(1))
    lg <- lg[slug %in% slugs]
    
    lg_ok  <- lg[status == "ok" & !is.na(timestamp_utc_parsed)]
    lg_pol <- lg_ok[policy_id == policy$policy_id]
    
    pick_latest <- function(dt) {
      dt[order(timestamp_utc_parsed)][, .SD[.N], by = slug]
    }
    
    latest <- if (nrow(lg_pol) > 0) pick_latest(lg_pol) else pick_latest(lg_ok)
    
    latest[, n_in  := as.integer(n_in)]
    latest[, n_out := as.integer(n_out)]
    latest[, dropped_total := n_in - n_out]
    latest[, kept_pct := round(100 * n_out / pmax(n_in, 1L), 1)]
    
    drop_cols <- grep("^dropped_", names(latest), value = TRUE)
    drop_cols_for_top <- setdiff(drop_cols, "dropped_total")
    
    if (length(drop_cols_for_top) > 0) {
      m <- as.matrix(latest[, ..drop_cols_for_top])
      suppressWarnings(storage.mode(m) <- "numeric")
      m[is.na(m)] <- 0
      
      max_n <- apply(m, 1, max)
      idx   <- max.col(m, ties.method = "first")
      
      latest[, top_drop_n := as.integer(max_n)]
      latest[, top_drop_reason := drop_cols_for_top[idx]]
      latest[top_drop_n == 0, top_drop_reason := NA_character_]
    } else {
      latest[, `:=`(top_drop_reason = NA_character_, top_drop_n = 0L)]
    }
    
    cat("\n============================================================\n")
    cat("Stage 03 summary (latest run per species)\n")
    cat("Policy preference:", policy$policy_id, "\n")
    cat("Runlog:", runlog_path, "\n")
    cat("============================================================\n\n")
    
    print(latest[, .(
      species, slug, policy_id,
      n_in, n_out, kept_pct, dropped_total,
      top_drop_reason, top_drop_n
    )][order(-dropped_total)])
    
    cat("\n-- Totals across these species (latest rows) --\n")
    print(latest[, .(
      n_in = sum(n_in, na.rm = TRUE),
      n_out = sum(n_out, na.rm = TRUE),
      dropped_total = sum(dropped_total, na.rm = TRUE),
      kept_pct = round(100 * sum(n_out, na.rm = TRUE) / pmax(sum(n_in, na.rm = TRUE), 1), 1)
    )])
    
    if (length(drop_cols) > 0) {
      reason_totals <- latest[, lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)), .SDcols = drop_cols]
      reason_totals <- melt(
        reason_totals,
        measure.vars = names(reason_totals),
        variable.name = "reason",
        value.name = "n_dropped"
      )[order(-n_dropped)]
      
      cat("\n-- Drop reasons (summed across these species; latest rows) --\n")
      print(reason_totals[n_dropped > 0][1:min(12, .N)])
    }
    
    cat("\nNote: if you encounter intermittent parsing issues when reading the runlog, it is safe to delete and regenerate it.\n")
    cat("This will not affect any processed datasets; it only recreates the log file.\n")
    cat("  # file.remove(runlog_path)\n")
  }
}
