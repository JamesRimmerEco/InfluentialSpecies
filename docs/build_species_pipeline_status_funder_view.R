# docs/build_species_pipeline_status_clear_view.R -------------------------------
#
# Purpose
#   Take the canonical internal pipeline status table (CSV) and produce a shorter,
#   commissioner-friendly “digest” CSV.
#
# Inputs
#   docs/derived/species_pipeline_status.csv
#
# Output
#   docs/derived/species_pipeline_status_clear_view.csv
#
# Optional output
#   docs/derived/species_pipeline_status_clear_view__glossary.txt
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Repo root ---------------------------------------------------------------
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

find_repo_root <- function(start_dir) {
  markers <- c(".git", "R", "data", "InfluentialSpecies.Rproj", "DESCRIPTION")
  d <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:15) {
    if (any(file.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root <- if (!is.null(this_file) && nzchar(this_file)) {
  find_repo_root(dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE)))
} else {
  find_repo_root(getwd())
}

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# Write a plain-text glossary alongside the CSV (overwritten on each run).
write_glossary <- TRUE

# ==============================================================================
# Input / output paths
# ==============================================================================

in_path  <- file.path(repo_root, "docs", "derived", "species_pipeline_status.csv")
out_path <- file.path(repo_root, "docs", "derived", "species_pipeline_status_clear_view.csv")

if (!file.exists(in_path)) {
  stop("Can't find canonical status table at: ", in_path, "\n",
       "Run docs/build_species_pipeline_status_csv.R first.")
}

# ---- Read canonical table ----------------------------------------------------
dt <- fread(in_path, na.strings = c("", "NA"))

# ---- Small helpers -----------------------------------------------------------
as_int <- function(x) suppressWarnings(as.integer(x))

pct1 <- function(num, den) {
  out <- rep(NA_real_, length(num))
  ok <- !is.na(num) & !is.na(den) & den > 0
  out[ok] <- round(100 * num[ok] / den[ok], 1)
  out
}

# ---- Pull + standardise expected fields -------------------------------------
dt[, n_gbif := as_int(n_gbif_after_strict_dedup)]
dt[, n_nbn  := as_int(n_nbn_after_strict_dedup)]
dt[, n_strict_total := as_int(n_total_after_strict_dedup)]
dt[, n_merge_final  := as_int(n_after_merge)]

dt[, n_filter_out := as_int(n_after_policy_filter)]

dt[, n_presence_cells := as_int(n_gridded_presence_cells)]
dt[, n_points_for_gee := as_int(n_gridded_points_for_gee)]

# ---- Derived metrics: attrition ---------------------------------------------
dt[, drop_merge_n := ifelse(!is.na(n_strict_total) & !is.na(n_merge_final),
                            n_strict_total - n_merge_final, NA_integer_)]
dt[, drop_merge_pct := pct1(drop_merge_n, n_strict_total)]

# Policy-drop needs a denominator; use n_before_policy_filter if present,
# but do not keep it in the output table.
dt[, n_filter_in := as_int(n_before_policy_filter)]
dt[, drop_policy_n := ifelse(!is.na(n_filter_in) & !is.na(n_filter_out),
                             n_filter_in - n_filter_out, NA_integer_)]
dt[, drop_policy_pct := pct1(drop_policy_n, n_filter_in)]

# ---- Progress placeholders (leave blank; can be filled later) ----------------
if (!("gee_assets_exist" %in% names(dt))) dt[, gee_assets_exist := ""]
if (!("gee_model_run" %in% names(dt))) dt[, gee_model_run := ""]

# ---- A readable stage summary -----------------------------------------------
dt[, pipeline_stage_reached := fifelse(
  !is.na(n_presence_cells), "Gridded",
  fifelse(!is.na(n_filter_out), "Filtered",
          fifelse(!is.na(n_merge_final), "Merged",
                  fifelse(!is.na(n_strict_total), "Downloaded", "Unknown")
          )
  )
)]

# ---- Final column selection + ordering --------------------------------------
out <- dt[, .(
  species_binomial,
  pipeline_stage_reached,
  
  n_gbif,
  n_nbn,
  n_strict_total,
  
  n_merge_final,
  drop_merge_n,
  drop_merge_pct,
  
  n_filter_out,
  drop_policy_n,
  drop_policy_pct,
  
  n_presence_cells,
  n_points_for_gee,
  
  gee_assets_exist,
  gee_model_run
)]

# ---- Sort: stage reached, then weakest (fewest cells) ------------------------
out[, stage_rank := fifelse(
  pipeline_stage_reached == "Gridded", 4L,
  fifelse(pipeline_stage_reached == "Filtered", 3L,
          fifelse(pipeline_stage_reached == "Merged", 2L,
                  fifelse(pipeline_stage_reached == "Downloaded", 1L, 0L)
          )
  )
)]

setorder(out, -stage_rank, n_presence_cells, na.last = TRUE)
out[, stage_rank := NULL]

# ---- Write output ------------------------------------------------------------
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
fwrite(out, out_path, na = "")

# ---- Glossary (overwrite each run) -------------------------------------------
if (isTRUE(write_glossary)) {
  gloss_path <- sub("\\.csv$", "__glossary.txt", out_path)
  
  gloss_lines <- c(
    "Species pipeline status — clear view (glossary)",
    "================================================",
    "",
    "This file defines each column in species_pipeline_status_clear_view.csv in plain terms.",
    "",
    "Columns",
    "-------",
    "",
    "species_binomial",
    "  Species scientific name (binomial) as listed in the project species list.",
    "",
    "pipeline_stage_reached",
    "  Furthest completed stage inferred from whether key outputs exist for that species:",
    "    Downloaded → Merged → Filtered → Gridded",
    "",
    "n_gbif",
    "  Number of GBIF occurrence records after GBIF-specific cleaning and strict de-duplication,",
    "  before any merging with NBN.",
    "",
    "n_nbn",
    "  Number of NBN occurrence records after NBN-specific cleaning and strict de-duplication,",
    "  before any merging with GBIF.",
    "",
    "n_strict_total",
    "  n_gbif + n_nbn: total occurrences from both sources before cross-source de-duplication.",
    "",
    "n_merge_final",
    "  Number of records after GBIF and NBN have been combined and duplicates removed between sources.",
    "",
    "drop_merge_n",
    "  Number of records removed during cross-source de-duplication:",
    "    drop_merge_n = n_strict_total - n_merge_final",
    "",
    "drop_merge_pct",
    "  Percentage of pre-merge records removed during cross-source de-duplication:",
    "    drop_merge_pct = (drop_merge_n / n_strict_total) * 100",
    "",
    "n_filter_out",
    "  Number of records remaining after applying the project’s policy filter.",
    "",
    "drop_policy_n",
    "  Number of records removed by policy filtering (difference between policy filter input and output).",
    "",
    "drop_policy_pct",
    "  Percentage removed by policy filtering:",
    "    drop_policy_pct = (drop_policy_n / policy_filter_input) * 100",
    "",
    "n_presence_cells",
    "  Number of grid cells containing at least one retained occurrence after gridding (unique occupied cells).",
    "",
    "n_points_for_gee",
    "  Number of point features written for modelling in Google Earth Engine.",
    "  Typically one point per occupied cell using the project’s first-observed-location-in-cell rule,",
    "  but can be lower if points are dropped by masking or missing predictor coverage.",
    "",
    "gee_assets_exist",
    "  Placeholder for tracking — whether the expected GEE assets exist for this species.",
    "",
    "gee_model_run",
    "  Placeholder for tracking — whether a model has been run for this species in GEE."
  )
  
  writeLines(gloss_lines, gloss_path)
}

# ---- Console summary ---------------------------------------------------------
cat("\n============================================================\n")
cat("Clear-view status table written\n")
cat("============================================================\n")
cat("Input:  ", in_path, "\n", sep = "")
cat("Output: ", out_path, "\n", sep = "")
if (isTRUE(write_glossary)) cat("Glossary: ", sub("\\.csv$", "__glossary.txt", out_path), "\n", sep = "")
cat("Species total: ", nrow(out), "\n", sep = "")
cat("Done.\n")
cat("============================================================\n\n")