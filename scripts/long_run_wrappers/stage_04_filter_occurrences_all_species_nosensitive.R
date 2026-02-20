# scripts/long_run_wrappers/stage_04_filter_occurrences_all_species_nosensitive.R
#
# Stage 04: Policy filtering for all species (pre-rasterisation)
#
# Purpose:
#   Run Stage 04 over the full authoritative InfluentialSpecies list using a single,
#   explicitly-defined policy and write per-species filtered outputs plus a runlog.
#
#   Stage 04 is policy, so every numerical decision lives here. This wrapper is a
#   long-run “home run” script: restart-safe, deterministic engine selection, and
#   designed to resume without losing progress.
#
# Sensitive species:
#   This run applies a single shared policy to all taxa (no sensitivity-specific overrides).
#   The Stage 04 engine supports per-species overrides via policy$per_species_overrides(),
#   but this wrapper does not use them.
#
# Inputs:
#   data/processed/03_qc_flagged/<slug>/occ_<slug>__qc_flagged.(parquet|rds)
#
# Outputs:
#   data/processed/04_filtered/<slug>/occ_<slug>__filtered.(parquet|rds)
#   data/processed/04_filtered/_runlog_04_filtered.csv
#
# ------------------------------------------------------------------------------

# ---- Find repo root (works when sourced from any scripts/ subfolder) ----------
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file) || !nzchar(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run via source('.../scripts/.../stage_04_filter_occurrences_all_species_nosensitive.R') from a file, not copy/paste."
  )
}

script_dir <- dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE))

find_repo_root <- function(start_dir) {
  markers <- c(".git", "R", "data", "InfluentialSpecies.Rproj", "DESCRIPTION")
  d <- start_dir
  for (i in 1:20) {
    if (any(file.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root <- find_repo_root(script_dir)

# ---- Load the Stage 04 engine (deterministic) --------------------------------
engine_fn <- file.path(repo_root, "R", "filter_occurrences.R")
if (!file.exists(engine_fn)) {
  stop("Can't find Stage 04 engine at: ", engine_fn)
}
source(engine_fn)

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# ---- Species list -------------------------------------------------------------
# Canonical Latin binomial list; one per row, first column.
species_csv <- file.path(repo_root, "data", "_meta", "species_list_binomial.csv")
if (!file.exists(species_csv)) stop("Can't find species list at: ", species_csv)

species_names <- read.csv(species_csv, stringsAsFactors = FALSE, header = TRUE)[[1]]
species_names <- as.character(species_names)
species_names <- trimws(species_names)
species_names <- species_names[!is.na(species_names) & nzchar(species_names)]
species_names <- species_names[tolower(species_names) != "binomial"]   # belt + braces
species_names <- unique(species_names)

if (length(species_names) < 10) {
  stop("Species list looks unexpectedly short (", length(species_names), "). Check: ", species_csv)
}

# ---- Stage 04 policy ----------------------------------------------------------
# Stage 04 is policy, so every numerical decision lives here.
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
  #
  # Sensitive species:
  #   No sensitivity-specific overrides are applied in this run.
  policy_id = "baseline_2000_unc1km_obs_plus_specimen_prov_gbif__nosensitive",
  
  # Keep only selected sources (NULL keeps all).
  keep_sources = NULL,              # e.g. c("GBIF", "NBN")
  
  # ---- Structural drops (usually TRUE) ----------------------------------------
  drop_missing_coords = TRUE,
  drop_coords_out_of_range = TRUE,
  drop_future_date = TRUE,
  
  # ---- Date completeness ------------------------------------------------------
  drop_missing_date = FALSE,
  
  # ---- Date window (applies to eventDate; may fallback to year if allowed) ----
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

in_root  <- file.path("data", "processed", "03_qc_flagged")
out_root <- file.path("data", "processed", "04_filtered")

# Long-run behaviour:
#   - overwrite=FALSE is restart-safe (skips species where output exists)
#   - set TRUE only when you intentionally want to rebuild Stage 04 outputs
overwrite     <- FALSE

write_parquet <- TRUE
write_rds     <- FALSE
write_runlog  <- TRUE

continue_on_error <- TRUE
verbose <- TRUE

# ==============================================================================
# RUN
# ==============================================================================

stage04_filter_occurrences(
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