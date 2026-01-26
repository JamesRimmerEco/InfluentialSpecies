# scripts/qc_flag_species_set_6sp_test.R ---------------------------------------
# Stage 02: QC flagging for a small species set
#
# Purpose:
#   Run the Stage 02 QC flagging step for a short species list, producing per-species
#   QC-flagged outputs and a run log.
#
# Outputs:
#   data/processed/02_qc_flagged/<slug>/occ_<slug>__qc_flagged.(parquet|rds)
#   data/processed/02_qc_flagged/_runlog_02_qc_flagged.csv
#
# Behaviour:
#   - Reads Stage 01 merged outputs from data/processed/01_merged.
#   - Adds QC flag columns and writes updated per-species files to Stage 02.
#   - Does not remove records; it only annotates them.
#
# Notes:
#   - Stage 01 inputs may be grouped or ungrouped. group_dir is used only for input discovery.
#   - Stage 02 outputs are always written ungrouped under data/processed/02_qc_flagged/<slug>/.
#   - Re-running later will rebuild only species whose Stage 01 inputs have changed if
#     refresh_if_inputs_newer=TRUE.

# ---- Find this script’s directory (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run this via source('.../scripts/qc_flag_species_set_6sp_test.R') from a file, not by copy/paste."
  )
}

script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

# ---- Load the workflow function ----
qc_fn <- file.path(repo_root, "R", "qc_flag_occurrences.R")
if (!file.exists(qc_fn)) {
  stop(
    "Can't find qc_flag_occurrences.R at: ", qc_fn,
    "\nCheck that the file exists at InfluentialSpecies/R/qc_flag_occurrences.R"
  )
}
source(qc_fn)

# ---- Key settings ------------------------------------------------------------
# If Stage 01 was written grouped, set e.g. "set_6sp_test".
# If Stage 01 was written ungrouped under data/processed/01_merged/<slug>/, set group_dir <- "".
group_dir <- ""

species_names <- c(
  "Myrmica sabuleti",
  "Myrmica scabrinodis",
  "Andrena fulva",
  "Sorex araneus",
  "Leptothorax acervorum",
  "Emberiza schoeniclus"
)

# ---- QC rule settings --------------------------------------------------------
# max_coord_uncertainty_m controls qc_flag_uncertainty_high.
max_coord_uncertainty_m <- 10000

# licence_expected is produced during the pull stage. If present and FALSE, it is flagged.
flag_if_unexpected_licence <- TRUE

# issues is mainly a GBIF concept. If present and non-empty, it is flagged.
flag_if_has_issues <- TRUE

# ---- Convenience: only run species where a Stage 01 merged input exists -------
only_run_if_stage01_exists <- TRUE

slugify_species <- function(species_name) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(species_name))
  slug <- gsub("^_+|_+$", "", slug)
  slug
}

normalise_group_dir <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) "" else x
}

find_stage01_merged_base <- function(repo_root, group_dir, slug) {
  group_dir <- normalise_group_dir(group_dir)
  base_name <- paste0("occ_", slug, "__merged")
  
  candidates <- character()
  
  if (nzchar(group_dir)) {
    candidates <- c(
      candidates,
      file.path(repo_root, "data", "processed", "01_merged", group_dir, slug, base_name)
    )
  }
  
  candidates <- c(
    candidates,
    file.path(repo_root, "data", "processed", "01_merged", slug, base_name)
  )
  
  candidates <- unique(candidates)
  
  for (b in candidates) {
    if (file.exists(paste0(b, ".parquet")) || file.exists(paste0(b, ".rds"))) return(b)
  }
  
  NA_character_
}

if (isTRUE(only_run_if_stage01_exists)) {
  keep <- vapply(species_names, function(sp) {
    slug <- slugify_species(sp)
    b <- find_stage01_merged_base(repo_root, group_dir, slug)
    !is.na(b)
  }, logical(1))
  
  dropped <- species_names[!keep]
  if (length(dropped) > 0) {
    message("Skipping species with no Stage 01 merged input present: ", paste(dropped, collapse = ", "))
  }
  species_names <- species_names[keep]
}

# ---- Run ---------------------------------------------------------------------
qc_flag_occurrences(
  species_names              = species_names,
  group_dir                  = group_dir,
  in_root                    = file.path("data", "processed", "01_merged"),
  out_root                   = file.path("data", "processed", "02_qc_flagged"),
  overwrite                  = TRUE, # Slower but necessary if upstream data are updated
  refresh_if_inputs_newer    = TRUE,
  continue_on_error          = TRUE,
  max_coord_uncertainty_m    = max_coord_uncertainty_m,
  flag_if_unexpected_licence = flag_if_unexpected_licence,
  flag_if_has_issues         = flag_if_has_issues,
  make_flag_count            = TRUE
)
