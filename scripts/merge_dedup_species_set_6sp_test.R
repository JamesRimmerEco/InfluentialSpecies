# scripts/merge_dedup_species_set_6sp_test.R ----------------------------------
# Stage 01: merge and safe 1-to-1 day duplicate auto-drop for a test species set
#
# Purpose:
#   Run the Stage 01 merge step for a small species list, producing per-species merged
#   outputs and a run log.
#
# Outputs:
#   data/processed/01_merged/<slug>/occ_<slug>__merged.(parquet|rds)
#   data/processed/01_merged/_runlog_01_merged.csv
#
# Behaviour:
#   - Reads raw-clean CSVs from data/raw/gbif and data/raw/nbn (created by pull_raw_occurrences()).
#   - Merges GBIF + NBN for each species and retains all original columns.
#   - Adds derived fields used downstream (IDs, parsed day, rounded coordinates, keys).
#   - Drops only the most conservative cross-source duplicates:
#       * exactly 1 GBIF + 1 NBN record share the same rounded coordinates and true day-level date
#       * the non-preferred source record is dropped (prefer_source)
#   - If only one source is available for a species, the merge still runs and writes output.
#
# Notes:
#   - Raw inputs may be grouped (data/raw/<src>/<group>/...) or ungrouped (data/raw/<src>/...).
#     This script can be used with either layout by setting group_dir appropriately.
#   - When GBIF downloads are still pending, some species may only have NBN inputs available.
#     Re-running later will incorporate new GBIF files if refresh_if_inputs_newer=TRUE.


# ---- Find this script’s directory (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run this via source('.../scripts/merge_dedup_species_set_6sp_test.R') from a file, not by copy/paste."
  )
}

script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

# ---- Load the workflow function ----
merge_fn <- file.path(repo_root, "R", "merge_dedup_occurrences.R")
if (!file.exists(merge_fn)) {
  stop(
    "Can't find merge_dedup_occurrences.R at: ", merge_fn,
    "\nCheck that the file exists at InfluentialSpecies/R/merge_dedup_occurrences.R"
  )
}
source(merge_fn)

# ---- Key settings ----
# If you used a grouped pull, set e.g. "set_6sp_test".
# If you wrote straight into data/raw/gbif and data/raw/nbn, set group_dir <- "".
group_dir <- ""   # <-- your current setup

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

# Optional: while GBIF downloads are pending, you can choose to only merge species
# that have at least one raw input present. (merge_occurrences() also handles missing
# inputs gracefully, so this is just a convenience.)
only_merge_if_any_input_exists <- TRUE

# ---- Helper: slugify + raw input discovery (matches merge script) ----
slugify_species <- function(species_name) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(species_name))
  slug <- gsub("^_+|_+$", "", slug)
  slug
}

normalise_group_dir <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) "" else x
}

find_raw_clean_csv <- function(repo_root, raw_dir, source, group_dir, slug) {
  src <- tolower(source)
  fname <- paste0(src, "_", slug, "_clean.csv")
  group_dir <- normalise_group_dir(group_dir)
  
  candidates <- c(
    file.path(repo_root, raw_dir, src, group_dir, slug, fname),
    file.path(repo_root, raw_dir, src, group_dir, fname),
    file.path(repo_root, raw_dir, src, slug, fname),
    file.path(repo_root, raw_dir, src, fname)
  )
  
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

if (isTRUE(only_merge_if_any_input_exists)) {
  raw_dir <- file.path("data", "raw")
  keep <- vapply(species_names, function(sp) {
    slug <- slugify_species(sp)
    gbif <- find_raw_clean_csv(repo_root, raw_dir, "gbif", group_dir, slug)
    nbn  <- find_raw_clean_csv(repo_root, raw_dir, "nbn",  group_dir, slug)
    (!is.na(gbif) && file.exists(gbif)) || (!is.na(nbn) && file.exists(nbn))
  }, logical(1))
  
  dropped <- species_names[!keep]
  if (length(dropped) > 0) {
    message("Skipping species with no raw inputs present: ", paste(dropped, collapse = ", "))
  }
  species_names <- species_names[keep]
}

# ---- Run ----
merge_occurrences(
  species_names           = species_names,
  group_dir               = group_dir,
  coord_round_dp          = coord_round_dp,
  prefer_source           = prefer_source,
  overwrite               = TRUE, # runs slower, but important if raw pull are updated (harder force than below argument)
  refresh_if_inputs_newer = TRUE,   # key: re-run later to pick up newly arrived GBIF downloads
  continue_on_error       = TRUE
)
