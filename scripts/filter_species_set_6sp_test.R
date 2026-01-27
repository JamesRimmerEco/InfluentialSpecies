# scripts/filter_species_set_6sp_test.R -----------------------------------------
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

# ---- Find this script’s directory (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run this via source('.../scripts/filter_species_set_6sp_test.R') from a file, not by copy/paste."
  )
}

script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

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

policy <- list(
  # Unique label for this exact set of filtering rules.
  # Change this whenever you change ANY threshold/switch (e.g. uncertainty 1 km -> 5 km),
  # so run logs and outputs can be traced back to the policy that produced them.
  policy_id = "baseline_2010_unc1km_in_situ_obs_only_gbif",
  
  # Keep only selected sources (NULL keeps all).
  keep_sources = NULL,              # e.g. c("GBIF", "NBN")
  
  # ---- Structural drops (usually TRUE) ----------------------------------------
  
  # Drop records with no usable coordinates (lon/lat missing).
  # These cannot be rasterised or mapped meaningfully.
  drop_missing_coords = TRUE,
  
  # Drop records with coordinates outside valid numeric ranges.
  # (Longitude must be between -180 and 180; latitude between -90 and 90.)
  drop_coords_out_of_range = TRUE,
  
  # Drop records dated in the future (event_day > today or year > current year).
  # These are almost always data entry/parse errors.
  drop_future_date = TRUE,
  
  # ---- Date completeness ------------------------------------------------------
  
  # TRUE drops rows with no event_day and no year.
  drop_missing_date = FALSE,
  
  # ---- Date window (applies to event_day; may fallback to year if allowed) ----
  
  min_date = "2000-01-01",          # set NA to disable
  max_date = NA,                    # set NA to disable
  
  # If event_day is missing but a year is available, allow records to pass
  # date filters using the year value only (e.g. year >= min_year).
  allow_year_only = TRUE,
  
  # If TRUE, drop any row without event_day regardless of year.
  require_event_day = FALSE,
  
  # Optional year window (only used when allow_year_only=TRUE).
  min_year = NA,                    # e.g. 2010
  max_year = NA,
  
  # ---- Coordinate uncertainty -------------------------------------------------
  
  # Maximum coordinate uncertainty (metres). Records with a reported uncertainty
  # greater than this are dropped. Set NA to disable uncertainty filtering.
  max_coord_uncertainty_m = 1000,   # 1 km; set 5000 for 5 km, 10000 for 10 km, etc.
  
  # Some records have no coordinateUncertaintyInMeters value (unknown uncertainty).
  # Choose whether to keep those records ("keep") or drop them as too uncertain ("drop").
  uncertainty_missing_action = "keep",
  
  # ---- GBIF issues handling ---------------------------------------------------
  
  # GBIF attaches an 'issues' string when its ingestion/interpretation pipeline flags
  # something about a record (e.g. COORDINATE_ROUNDED, COUNTRY_DERIVED_FROM_COORDINATES).
  # Many issue codes are common and harmless, so the default is to ignore here as we pick up
  # issues ourselves anyway.
  #
  # Options:
  #   "ignore"         : keep all records regardless of issues (recommended default)
  #   "drop_any"       : drop GBIF records whenever issues is non-empty (very strict)
  #   "drop_blacklist" : drop GBIF records only if issues contains one or more specified codes
  gbif_issues_mode = "ignore",
  
  # Only used when gbif_issues_mode = "drop_blacklist".
  # Provide exact GBIF issue codes to exclude. Leave empty to exclude none. Remember that many
  # of these issues will be picked up ourselves anyway.
  issues_blacklist = c(
    # Examples:
    # "ZERO_COORDINATE",
    # "COORDINATE_OUT_OF_RANGE",
    # "RECORDED_DATE_INVALID"
  ),
  
  # ---- Licence handling -------------------------------------------------------
  
  # Only applies if licence_expected exists. TRUE drops records where licence_expected == FALSE.
  drop_unexpected_licence = FALSE,
  
  # ---- basisOfRecord handling (optional; applies if basisOfRecord exists) -----
  #
  # NOTE:
  # NBN often has basisOfRecord == NA, so we do NOT use allowed_basis_of_record/drop_basis_of_record
  # globally here (that would accidentally drop large amounts of NBN).
  # Instead we enforce an explicit "in situ observations only" rule for GBIF in extra_drop_rules below.
  
  allowed_basis_of_record = NULL,   # e.g. c("HUMAN_OBSERVATION","OBSERVATION","MACHINE_OBSERVATION")
  drop_basis_of_record = NULL,      # e.g. c("FOSSIL_SPECIMEN","PRESERVED_SPECIMEN")
  
  # ---- Taxon rank gating (optional; applies if taxonRank exists) --------------
  
  allowed_taxon_rank = NULL,        # e.g. c("SPECIES")
  
  # ---- NBN certainty gating (optional; applies if source=="NBN" and column exists)
  
  # Column name must match your merged schema (e.g. "certainty" or similar).
  nbn_certainty_col = NULL,         # e.g. "certainty"
  nbn_allowed_certainty = NULL,     # e.g. c("definite","probable")
  
  # ---- Extra rules (advanced) -------------------------------------------------
  
  # Named list of functions(dt) -> logical drop vector (TRUE means drop).
  extra_drop_rules = list(
    # In situ / contemporary: drop any GBIF record that is not clearly an observation.
    #
    # We KEEP only:
    #   - HUMAN_OBSERVATION
    #   - OBSERVATION
    #   - MACHINE_OBSERVATION
    #
    # We DROP anything else (e.g. PRESERVED_SPECIMEN, FOSSIL_SPECIMEN, LIVING_SPECIMEN,
    # MATERIAL_SAMPLE, MATERIAL_CITATION, OCCURRENCE, etc.)
    drop_gbif_non_observation_basis = function(dt) {
      if (!("source" %in% names(dt)) || !("basisOfRecord" %in% names(dt))) {
        return(rep(FALSE, nrow(dt)))
      }
      
      bor <- toupper(trimws(as.character(dt$basisOfRecord)))
      allowed <- c("HUMAN_OBSERVATION", "OBSERVATION", "MACHINE_OBSERVATION")
      
      dt$source == "GBIF" & !(bor %in% allowed)
    }
    
    # Example:
    # drop_pre2000_gbif = function(dt) {
    #   yr <- suppressWarnings(as.integer(dt$year))
    #   dt$source == "GBIF" & !is.na(yr) & yr < 2000
    # }
  )
)

# ==============================================================================
# RUN SETTINGS
# ==============================================================================

in_root  <- file.path("data", "processed", "02_qc_flagged")
out_root <- file.path("data", "processed", "03_filtered")

# Stage 03 is expected to be re-run often as policy changes.
overwrite     <- TRUE
write_parquet <- TRUE
write_rds     <- FALSE

# Optional: keep a compact run history (recommended for PI "what changed?" questions)
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
