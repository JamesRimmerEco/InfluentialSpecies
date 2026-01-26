# InfluentialSpecies/scripts/pull_raw_species_set_6sp_test_v_2.0.R
#
# Run the raw occurrence pull for a defined set of species.
# Folder paths + file saving are handled inside: InfluentialSpecies/R/pull_raw_occurrences.R
#
# Outputs write into GBIF/NBN subfolders automatically, e.g.:
#   data/raw/gbif/<species_slug>/gbif_<species_slug>_clean.csv
#   data/raw/nbn/<species_slug>/nbn_<species_slug>_clean.csv
# 

# ---- Find this script’s directory (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop("Can't determine script path (sys.frame(1)$ofile is NULL). ",
       "Run this via source('.../scripts/pull_raw_species_set_6sp_test.R') from a file, not by copy/paste.")
}

script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

# ---- Load the workflow function ----
pull_fn <- file.path(repo_root, "R", "pull_raw_occurrences.R")
if (!file.exists(pull_fn)) {
  stop("Can't find pull_raw_occurrences.R at: ", pull_fn,
       "\nCheck that the file exists at InfluentialSpecies/R/pull_raw_occurrences.R")
}
source(pull_fn)

# ---- Key settings ----
nbn_email <- "jamesrimmer92@mail.com"  # email used for NBN (galah) access

# If we previously ran this before the QA columns were added, we keep TRUE:
# the pull script will auto-repull if cached CSVs are missing the new QA columns.
#
# Patch note:
#   Stage 00 now also retains record-type / provenance fields (e.g. basisOfRecord, taxonRank,
#   occurrenceStatus, datasetKey, etc.). Older cached CSVs missing these will be treated as
#   stale automatically and re-pulled when the upstream APIs are available.
use_cache <- TRUE

# write straight into data/raw/gbif and data/raw/nbn
group_dir <- ""  

# Species list (Latin names only; common names in comments for readability but not needed)
species_names <- c(
  "Myrmica sabuleti",        # Myrmica sabuleti
  "Myrmica scabrinodis",     # Scabrous ant
  "Andrena fulva",           # Tawny mining bee
  "Sorex araneus",           # Common shrew
  "Leptothorax acervorum",   # Small red ant
  "Emberiza schoeniclus"     # Reed bunting
)

# ---- Run ----
pull_raw_occurrences(
  species_names      = species_names,
  group_dir          = group_dir,
  nbn_email          = nbn_email,
  use_cache          = use_cache,
  species_subdir     = FALSE,
  gbif_method        = "auto",
  gbif_download_wait = FALSE
)
