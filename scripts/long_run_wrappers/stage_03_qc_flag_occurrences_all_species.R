# scripts/long_run_wrappers/stage_03_qc_flag_occurrences_all_species.R
#
# Stage 03: QC flagging for all species (resume safely)
#
# Purpose:
#   Run Stage 03 QC flagging over the full authoritative InfluentialSpecies list.
#   This stage annotates records with QC flags but does not drop records.
#
# Inputs:
#   data/processed/02_merged/<slug>/occ_<slug>__merged.(parquet|rds)
#   (also supports grouped layout: data/processed/02_merged/<group>/<slug>/...)
#
# Outputs:
#   data/processed/03_qc_flagged/<slug>/occ_<slug>__qc_flagged.(parquet|rds)
#   data/processed/03_qc_flagged/_runlog_03_qc_flagged.csv
#
# ------------------------------------------------------------------------------

# ---- Find repo root (works when sourced from any scripts/ subfolder) ----------

this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (is.null(this_file) || !nzchar(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run via source('.../scripts/long_run_wrappers/stage_03_qc_flag_occurrences_all_species.R') from a file, not copy/paste."
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

# Engines in this repo locate the project via getwd(); make this deterministic.
setwd(repo_root)

# ---- Load Stage 03 engine (deterministic) ------------------------------------

qc_candidates <- c(
  file.path(repo_root, "R", "qc_flag_occurrences.R")
)

qc_fn <- qc_candidates[file.exists(qc_candidates)][1]

if (is.na(qc_fn) || !nzchar(qc_fn)) {
  stop(
    "Can't find Stage 03 engine. Checked:\n  - ",
    paste(qc_candidates, collapse = "\n  - ")
  )
}

source(qc_fn)

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# ---- Species list -------------------------------------------------------------

# Canonical Latin binomial list; one per row, first column (no header).
species_csv <- file.path(repo_root, "data", "_meta", "species_list_binomial.csv")
if (!file.exists(species_csv)) stop("Can't find species list at: ", species_csv)

species_names <- read.csv(species_csv, stringsAsFactors = FALSE, header = FALSE)[[1]]
species_names <- as.character(species_names)
species_names <- trimws(species_names)
species_names <- species_names[!is.na(species_names) & nzchar(species_names)]
species_names <- species_names[tolower(species_names) != "binomial"]   # belt + braces
species_names <- unique(species_names)

if (length(species_names) < 10) {
  stop("Species list looks unexpectedly short (", length(species_names), "). Check: ", species_csv)
}

# ==============================================================================
# RUN SETTINGS
# ==============================================================================

# Stage 02 inputs (merged)
in_root <- file.path("data", "processed", "02_merged")

# Stage 03 outputs (qc flagged)
out_root <- file.path("data", "processed", "03_qc_flagged")

# If your Stage 02 outputs are under a group directory (e.g. data/processed/02_merged/<group>/<slug>/...),
# set group_dir to that folder name; otherwise leave as "" for ungrouped.
group_dir <- ""

# Long-run behaviour:
#   - overwrite=FALSE is restart-safe (skips species where output exists)
#   - refresh_if_inputs_newer=TRUE rebuilds outputs only when Stage 02 inputs are newer
overwrite <- FALSE
refresh_if_inputs_newer <- TRUE
continue_on_error <- TRUE

# QC flag parameters (Stage 03 does not drop records)
max_coord_uncertainty_m <- 10000
flag_if_unexpected_licence <- TRUE
flag_if_has_issues <- TRUE
make_flag_count <- TRUE

# ==============================================================================
# RUN
# ==============================================================================

qc_flag_occurrences(
  species_names = species_names,
  group_dir = group_dir,
  in_root = in_root,
  out_root = out_root,
  overwrite = overwrite,
  refresh_if_inputs_newer = refresh_if_inputs_newer,
  continue_on_error = continue_on_error,
  max_coord_uncertainty_m = max_coord_uncertainty_m,
  flag_if_unexpected_licence = flag_if_unexpected_licence,
  flag_if_has_issues = flag_if_has_issues,
  make_flag_count = make_flag_count
)