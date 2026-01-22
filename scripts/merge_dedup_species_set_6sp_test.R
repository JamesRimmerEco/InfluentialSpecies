# scripts/merge_dedup_species_set_6sp_test.R ----------------------------------
# Stage 01: merge and safe 1-to-1 day duplicate auto-drop for 6 species test set
#
# Notes:
# - Uses InfluentialSpecies/R/merge_dedup_occurrences.R
# - Input discovery is robust (grouped/ungrouped, subdir/flat) in the merge function.
# - If GBIF downloads are pending for a species, you can still run this:
#     * it will merge whatever is available now (often NBN-only),
#     * then on later runs it will automatically refresh when new raw CSVs appear,
#       as long as refresh_if_inputs_newer = TRUE (default here).

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
# If you wrote straight into data/raw/gbif and data/raw/nbn, set "".
group_dir <- "set_6sp_test"
# group_dir <- ""   # <-- quick switch for ad-hoc / ungrouped test pulls

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

# Optional: while GBIF downloads are pending, you may want to merge only
# species that currently have at least one raw input present (GBIF or NBN).
# Set TRUE to skip species with no current inputs at all.
only_merge_if_any_input_exists <- TRUE

# ---- Helper: quick existence check for raw inputs (supports old/new layouts) ----
# (merge_occurrences() will also handle missing inputs gracefully; this just helps
# you avoid cluttering logs with "no_inputs" while a species is pending.)
normalise_group_dir <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) "" else x
}

slugify_species <- function(species_name) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(species_name))
  slug <- gsub("^_+|_+$", "", slug)
  slug
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
    any(file.exists(c(gbif, nbn)))
  }, logical(1))
  
  if (!any(keep)) {
    message("No raw inputs found yet for any species in this set. Nothing to merge.")
    quit(save = "no", status = 0)
  }
  
  skipped <- species_names[!keep]
  if (length(skipped) > 0) {
    message("Skipping species with no current inputs (yet):")
    for (s in skipped) message(" - ", s)
  }
  
  species_names <- species_names[keep]
}

# ---- Run ----
merge_occurrences(
  species_names           = species_names,
  group_dir               = group_dir,
  coord_round_dp          = coord_round_dp,
  prefer_source           = prefer_source,
  overwrite               = FALSE,
  refresh_if_inputs_newer = TRUE,
  continue_on_error       = TRUE
)
