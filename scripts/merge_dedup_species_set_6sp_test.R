# scripts/merge_dedup_species_set_6sp_test.R ----------------------------------
# Stage 01: merge and safe 1-to-1 day duplicate auto-drop for 6 species test set

# ---- Find this script’s directory (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop("Can't determine script path (sys.frame(1)$ofile is NULL). ",
       "Run this via source('.../scripts/merge_dedup_species_set_6sp_test.R') from a file, not by copy/paste.")
}

script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

# ---- Load the workflow function ----
merge_fn <- file.path(repo_root, "R", "merge_dedup_occurrences.R")
if (!file.exists(merge_fn)) {
  stop("Can't find merge_dedup_occurrences.R at: ", merge_fn,
       "\nCheck that the file exists at InfluentialSpecies/R/merge_dedup_occurrences.R")
}
source(merge_fn)

# ---- Key settings ----
group_dir <- "set_6sp_test"

species_names <- c(
  "Myrmica sabuleti",
  "Myrmica scabrinodis",
  "Andrena fulva",
  "Sorex araneus",
  "Leptothorax acervorum",
  "Emberiza schoeniclus"
)

# Conservative duplicate rule settings
coord_round_dp <- 4
prefer_source  <- "GBIF"

# ---- Run ----
merge_occurrences(
  species_names     = species_names,
  group_dir         = group_dir,
  species_subdir    = TRUE,
  coord_round_dp    = coord_round_dp,
  prefer_source     = prefer_source,
  overwrite         = FALSE,
  continue_on_error = TRUE
)
