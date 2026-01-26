# InfluentialSpecies — pull raw occurrences (GBIF + NBN Atlas) -------------------
#
# Purpose:
#   Pull occurrence data for one or more species from:
#     - GBIF (Europe-wide scope, via rgbif)
#     - NBN Atlas (UK-only, via galah)
#
#   Apply only very basic screening suitable for raw outputs:
#     - keep only records with usable coordinates (no geographical context is not useful)
#     - light-touch within-source de-duplication only (exact duplicates):
#         GBIF: duplicate gbifID, and exact repeats of (lon, lat, date)
#         NBN : duplicate recordID, and exact repeats of (lon, lat, date)
#
#    Save outputs for each species:
#    If species_subdir = TRUE: data/raw/gbif/<group_dir>/<slug>/gbif_<slug>_clean.csv
#
#    If species_subdir = FALSE: data/raw/gbif/<group_dir>/gbif_<slug>_clean.csv (or directly data/raw #    /gbif/... if group_dir is blank)
#
#   Save checkpoints to:
#     data/_checkpoints/gbif/gbif_pull_checkpoint_<slug>.rds
#     data/_checkpoints/nbn/nbn_pull_checkpoint_<slug>.rds
#
# How GBIF pulls work (important):
#   - GBIF "search" (occ_search) is hard-limited to 100,000 records per query.
#   - This script checks the expected GBIF record count per species:
#       * If <= 100,000: it uses occ_search paging and writes the CSV immediately.
#       * If > 100,000: it uses GBIF downloads (occ_download), which are asynchronous.
#   - For >100k species, the first run typically:
#       * submits a download job to GBIF,
#       * saves the download key in the checkpoint,
#       * skips to the next species (so the whole run doesn't stall),
#       * prints a final warning listing any species still pending.
#   - When you run the same pull script again later, it automatically:
#       * reads the saved download key from the checkpoint,
#       * checks whether the download has finished,
#       * fetches/unzips/cleans the data once ready, and writes the final CSV.
#   - Once a species is complete (CSV exists AND checkpoint is marked complete),
#     re-running does NOT re-download or re-pull that species; it is treated as cached.
#
# Typical usage:
#   - You usually run this via a wrapper script in /scripts/ (e.g. pull_raw_species_set_*.R)
#   - If no species exceed 100k, one run is enough.
#   - If any species exceed 100k, you may need to re-run the wrapper script one or more
#     times until the final "GBIF WARNING" list disappears (all downloads completed).
#
# QA / certainty fields included in the "raw clean" outputs (NEW):
#   GBIF:
#     - coordinateUncertaintyInMeters
#     - identificationVerificationStatus
#     - issues
#     - identifiedBy
#     - dateIdentified
#   NBN:
#     - coordinateUncertaintyInMeters
#     - identificationVerificationStatus
#     - identifiedBy
#     - coordinatePrecision (often NA, but included if selectable)
#
# Licence handling:
#   We allow all licences (no licence-based filtering).
#   However, we still define "expected" licence sets for each source, and if any
#   other licence types appear we:
#     (i) flag this clearly to the console, and
#     (ii) write a per-species log file under data/raw/licence_flags/
#   These log files are only created when unexpected licence types are
#   present (if there are none, this will never appear).
#
# Notes:
#   - If cached *_clean.csv files were created before the QA fields were added,
#     this script will automatically treat them as stale and re-pull (probably only
#     will occur for wasps, or if we decide to add more columns).
# ------------------------------------------------------------------------------ 


suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(lubridate)
  library(tibble)
  library(rgbif)
  library(galah)
})

# ---- Helper: load local (gitignored) credentials if present -------------------
load_local_gbif_credentials <- function(repo_root) {
  candidates <- c(
    file.path(repo_root, "credentials.R"),
    file.path(repo_root, "data", "credentials.R")
  )
  
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) return(invisible(FALSE))
  cred_file <- existing[1]
  
  cred_env <- new.env(parent = baseenv())
  sys.source(cred_file, envir = cred_env)
  
  if (exists("GBIF_USER", envir = cred_env, inherits = FALSE) &&
      exists("GBIF_PWD",  envir = cred_env, inherits = FALSE) &&
      exists("GBIF_EMAIL",envir = cred_env, inherits = FALSE)) {
    
    Sys.setenv(
      GBIF_USER  = get("GBIF_USER",  envir = cred_env, inherits = FALSE),
      GBIF_PWD   = get("GBIF_PWD",   envir = cred_env, inherits = FALSE),
      GBIF_EMAIL = get("GBIF_EMAIL", envir = cred_env, inherits = FALSE)
    )
    return(invisible(TRUE))
  }
  
  warning("Found credentials file but it did not define GBIF_USER / GBIF_PWD / GBIF_EMAIL.")
  invisible(FALSE)
}

# ---- Helper: find repo root robustly -----------------------------------------
get_repo_root <- function() {
  wd <- getwd()
  
  if (dir.exists(file.path(wd, "data"))) return(wd)
  if (dir.exists(file.path(wd, "..", "data"))) return(normalizePath(file.path(wd, ".."), mustWork = FALSE))
  
  stop(
    "Can't locate repo root.\n",
    "Expected to find a 'data/' folder at either:\n",
    "  - ", file.path(wd, "data"), "\n",
    "  - ", file.path(wd, "..", "data"), "\n",
    "Set your working directory to the project root (InfluentialSpecies) and try again."
  )
}

repo_root <- get_repo_root()
load_local_gbif_credentials(repo_root)

# ---- Helper: slugify a species name ------------------------------------------
slugify_species <- function(species_name) {
  slug <- str_replace_all(tolower(species_name), "[^a-z0-9]+", "_")
  slug <- str_replace_all(slug, "^_+|_+$", "")
  slug
}

# "Slugify" - just means creating a user-friendly descriptive domain, easier to read
# than codes.

# ---- Helper: Make group_dir robust when blank --------------------------------
normalise_group_dir <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) "" else x
}

# ---- Helper: write unexpected licence log (only if needed) -------------------
write_unexpected_licence_log <- function(species_name, slug, source_name, unexpected_tbl, repo_root) {
  if (nrow(unexpected_tbl) == 0) return(invisible(NULL))
  
  log_dir <- file.path(repo_root, "data", "raw", "licence_flags")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  log_file <- file.path(log_dir, paste0(tolower(source_name), "_", slug, "_unexpected_licences.csv"))
  readr::write_csv(unexpected_tbl, log_file)
  
  message(
    "\n[LICENCE FLAG] Unexpected licence types detected for ",
    source_name, " (", species_name, ").\n",
    "  Wrote log: ", log_file
  )
  
  invisible(log_file)
}

# ==============================================================================
# GBIF pull (Europe-wide) ------------------------------------------------------
# ==============================================================================

pull_gbif_clean <- function(species_name,
                            region_scope = "EUROPE",
                            group_dir = "",
                            species_subdir = FALSE,
                            pause_s = 0.25,
                            page_size = 1000,
                            max_records = Inf,
                            expected_licences_gbif = c("CC0_1_0", "CC_BY_4_0", "CC_BY_NC_4_0"),
                            use_cache = TRUE,
                            gbif_method = c("auto", "search", "download"),
                            gbif_download_wait = FALSE,
                            gbif_search_hard_limit = 100000L,
                            gbif_user = Sys.getenv("GBIF_USER"),
                            gbif_pwd = Sys.getenv("GBIF_PWD"),
                            gbif_email = Sys.getenv("GBIF_EMAIL")) {
  
  repo_root <- get_repo_root()
  group_dir <- normalise_group_dir(group_dir)
  gbif_out_root <- if (nzchar(group_dir)) {
    file.path(repo_root, "data", "raw", "gbif", group_dir)
  } else {
    file.path(repo_root, "data", "raw", "gbif")
  }
  
  slug <- slugify_species(species_name)
  
  # Output dirs
  ckpt_root     <- file.path(repo_root, "data", "_checkpoints")
  gbif_ckpt_dir <- file.path(ckpt_root, "gbif")
  
  dir.create(gbif_out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(gbif_ckpt_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Species subfolder
  gbif_out_dir <- if (isTRUE(species_subdir)) file.path(gbif_out_root, slug) else gbif_out_root
  dir.create(gbif_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Clean output file (used as a cache when re-running)
  gbif_outfile <- file.path(gbif_out_dir, paste0("gbif_", slug, "_clean.csv"))
  
  # Checkpoint file (also stores download keys + completion state)
  ckpt_file <- file.path(gbif_ckpt_dir, paste0("gbif_pull_checkpoint_", slug, ".rds"))
  
  # Licence normaliser (GBIF uses URLs with 'licenses' in the path, but can accept codes as well)
  lic_normalise_gbif <- function(x) {
    x_chr <- as.character(x)
    x_up  <- toupper(x_chr)
    
    dplyr::case_when(
      is.na(x_chr) ~ NA_character_,
      
      # already-normalised codes (seen in download files)
      x_up %in% c("CC0_1_0", "CC_BY_4_0", "CC_BY_NC_4_0") ~ x_up,
      
      # URL forms (seen in occ_search)
      stringr::str_detect(x_chr, "publicdomain/zero/1.0") ~ "CC0_1_0",
      stringr::str_detect(x_chr, "licenses/by-nc/4.0") ~ "CC_BY_NC_4_0",
      stringr::str_detect(x_chr, "licenses/by/4.0") ~ "CC_BY_4_0",
      
      TRUE ~ NA_character_
    )
  }
  
  
  # Columns that MUST exist in cached outputs to be considered "current schema"
  required_cache_cols <- c(
    "coordinateUncertaintyInMeters",
    "identificationVerificationStatus",
    "issues",
    # Record type / provenance fields (used for Stage 03+ filtering/diagnostics)
    "basisOfRecord",
    "taxonRank",
    "occurrenceStatus",
    "datasetKey",
    "datasetName",
    "publishingOrgKey",
    "institutionCode",
    "collectionCode"
  )
  
  # Always return a tibble with the expected schema (even if 0 rows)
  empty_gbif_clean <- function() {
    tibble::tibble(
      source = character(),
      species = character(),
      gbifID = character(),
      occurrenceID = character(),
      lon = numeric(),
      lat = numeric(),
      date = character(),
      year = integer(),
      country = character(),
      licence_raw = character(),
      licence = character(),
      licence_expected = logical(),
      coordinateUncertaintyInMeters = numeric(),
      identificationVerificationStatus = character(),
      issues = character(),
      identifiedBy = character(),
      dateIdentified = character(),
      # Record type / provenance fields
      basisOfRecord = character(),
      taxonRank = character(),
      occurrenceStatus = character(),
      datasetKey = character(),
      datasetName = character(),
      publishingOrgKey = character(),
      institutionCode = character(),
      collectionCode = character()
    )
  }
  
  # Load / initialise checkpoint (backward compatible with older ckpt schema)
  ckpt <- list(
    schema_version = 2,
    mode = NA_character_,          # "search" or "download"
    complete = FALSE,              # TRUE only once we have a full dataset saved
    total_expected = NA_integer_,
    # search paging state
    start = 0,
    all_pages = list(),
    # download state
    download_key = NA_character_,
    download_status = NA_character_,
    last_updated = as.character(Sys.time())
  )
  
  if (file.exists(ckpt_file)) {
    old <- readRDS(ckpt_file)
    ckpt <- utils::modifyList(ckpt, old)
  }
  
  # ---------------------------------------------------------------------------
  # GBIF taxon resolution
  # ---------------------------------------------------------------------------
  bb <- name_backbone(name = species_name)
  taxon_key <- bb$usageKey
  message("GBIF match: ", bb$scientificName, " (usageKey=", taxon_key, ", matchType=", bb$matchType, ")")
  
  # ---------------------------------------------------------------------------
  # Cache check (only trusted as "complete" if checkpoint says complete=TRUE)
  # ---------------------------------------------------------------------------
  use_cached <- isTRUE(use_cache) && file.exists(gbif_outfile)
  if (use_cached) {
    gbif_cached <- readr::read_csv(gbif_outfile, show_col_types = FALSE)
    missing_cols <- setdiff(required_cache_cols, names(gbif_cached))
    
    if (length(missing_cols) > 0) {
      message(
        "Found existing GBIF clean file but it is missing required columns for the current schema: ",
        paste(missing_cols, collapse = ", "),
        "\nRe-pulling from GBIF: ", gbif_outfile
      )
      use_cached <- FALSE
    } else if (isTRUE(ckpt$complete)) {
      message("Found existing GBIF clean file + checkpoint marked complete; skipping API pull: ", gbif_outfile)
      gbif_clean <- gbif_cached %>%
        mutate(
          licence = as.character(licence),
          licence_raw = as.character(licence_raw),
          licence_expected = !is.na(licence) & licence %in% expected_licences_gbif,
          identificationVerificationStatus = as.character(identificationVerificationStatus),
          issues = as.character(issues),
          basisOfRecord = as.character(basisOfRecord),
          taxonRank = as.character(taxonRank),
          occurrenceStatus = as.character(occurrenceStatus),
          datasetKey = as.character(datasetKey),
          datasetName = as.character(datasetName),
          publishingOrgKey = as.character(publishingOrgKey),
          institutionCode = as.character(institutionCode),
          collectionCode = as.character(collectionCode)
        )
      attr(gbif_clean, "gbif_status") <- list(
        state = "complete",
        method = if (!is.null(ckpt$mode) && nzchar(ckpt$mode)) ckpt$mode else "unknown"
      )
      return(gbif_clean)
    } else {
      message("Found existing GBIF clean file, but checkpoint is not marked complete; verifying completeness...")
      # We'll fall through and decide what to do based on record count + method
      use_cached <- FALSE
    }
  }
  
  # ---------------------------------------------------------------------------
  # Decide method: search (<=100k) vs download (>100k)
  #   - occ_search is capped at 100,000 records per query, so >100k requires downloads.
  # ---------------------------------------------------------------------------
  gbif_method <- match.arg(gbif_method)
  
  # Determine expected count (cheap; does not pull the data)
  total_expected <- tryCatch(
    rgbif::occ_count(
      taxonKey = taxon_key,
      continent = region_scope,
      hasCoordinate = TRUE
    ),
    error = function(e) NA_integer_
  )
  
  if (is.na(total_expected)) {
    message("[GBIF] Warning: could not determine expected count (occ_count failed). Defaulting to method=", gbif_method)
  } else {
    message("[GBIF] Expected rows (query): ", total_expected)
  }
  
  if (gbif_method == "auto") {
    if (!is.na(total_expected) && total_expected > gbif_search_hard_limit) {
      gbif_method <- "download"
    } else {
      gbif_method <- "search"
    }
  }
  
  # If user forced "search" but count suggests it will be capped, switch to download (to ensure full data).
  if (gbif_method == "search" && !is.na(total_expected) && total_expected > gbif_search_hard_limit) {
    message("[GBIF] Count exceeds ", gbif_search_hard_limit, ". occ_search() cannot retrieve >100k; switching to downloads.")
    gbif_method <- "download"
  }
  
  # Helper for safe credential check
  have_gbif_creds <- function() {
    nzchar(gbif_user) && nzchar(gbif_pwd) && nzchar(gbif_email)
  }
  
  # ---------------------------------------------------------------------------
  # DOWNLOAD path (unlimited, async) — submit now, resume on next run
  # ---------------------------------------------------------------------------
  if (gbif_method == "download") {
    
    ckpt$mode <- "download"
    ckpt$total_expected <- total_expected
    ckpt$last_updated <- as.character(Sys.time())
    saveRDS(ckpt, ckpt_file)
    
    if (!have_gbif_creds()) {
      message(
        "\n[GBIF][INCOMPLETE] This species requires GBIF downloads (>100k expected), ",
        "but GBIF credentials are not available in environment variables.\n",
        "Set GBIF_USER, GBIF_PWD, GBIF_EMAIL (e.g., in ~/.Renviron), then re-run.\n"
      )
      gbif_clean <- empty_gbif_clean()
      attr(gbif_clean, "gbif_status") <- list(state = "needs_credentials", method = "download", expected = total_expected)
      return(gbif_clean)
    }
    
    # If no key yet, submit a download and return (so the run can continue to other species)
    if (is.null(ckpt$download_key) || is.na(ckpt$download_key) || !nzchar(ckpt$download_key)) {
      
      dl <- rgbif::occ_download(
        rgbif::pred_and(
          rgbif::pred("taxonKey", taxon_key),
          rgbif::pred("continent", region_scope),
          rgbif::pred("hasCoordinate", TRUE)
        ),
        user = gbif_user, pwd = gbif_pwd, email = gbif_email
      )
      
      # occ_download() may return either a list with $key OR an atomic key (character)
      dl_key <- if (is.list(dl) && "key" %in% names(dl)) dl$key else as.character(dl)
      
      ckpt$download_key <- dl_key
      ckpt$download_status <- "SUBMITTED"
      ckpt$complete <- FALSE
      ckpt$last_updated <- as.character(Sys.time())
      saveRDS(ckpt, ckpt_file)
      
      message(
        "\n[GBIF] Download submitted for ", species_name, ".\n",
        "  key: ", ckpt$download_key, "\n",
        "  Status: SUBMITTED\n",
        "Re-run the script later to resume and fetch the finished download.\n"
      )
      
      gbif_clean <- empty_gbif_clean()
      attr(gbif_clean, "gbif_status") <- list(state = "pending_download", method = "download", key = ckpt$download_key, expected = total_expected)
      return(gbif_clean)
    }
    
    # We have a download key — check status, optionally wait, then fetch
    key <- ckpt$download_key
    
    if (isTRUE(gbif_download_wait)) {
      message("[GBIF] Waiting for download to finish (key=", key, ") ...")
      rgbif::occ_download_wait(key, status_ping = 30, quiet = FALSE)
    }
    
    meta <- tryCatch(rgbif::occ_download_meta(key), error = function(e) NULL)
    status <- if (!is.null(meta) && !is.null(meta$status)) meta$status else NA_character_
    ckpt$download_status <- status
    ckpt$last_updated <- as.character(Sys.time())
    saveRDS(ckpt, ckpt_file)
    
    if (is.na(status) || !identical(status, "SUCCEEDED")) {
      message(
        "\n[GBIF] Download not ready yet for ", species_name, " (key=", key, ", status=", status, ").\n",
        "Skipping for now — re-run later to resume.\n"
      )
      gbif_clean <- empty_gbif_clean()
      attr(gbif_clean, "gbif_status") <- list(state = "pending_download", method = "download", key = key, status = status, expected = total_expected)
      return(gbif_clean)
    }
    
    # Fetch zip
    zip_path <- rgbif::occ_download_get(key, path = gbif_ckpt_dir, overwrite = TRUE)
    if (is.list(zip_path) && "path" %in% names(zip_path)) zip_path <- zip_path$path
    
    tmpdir <- tempfile("gbif_dwc_")
    dir.create(tmpdir)
    utils::unzip(zip_path, exdir = tmpdir)
    
    occ_file <- list.files(
      tmpdir,
      pattern = "occurrence\\.(txt|csv)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )[1]
    
    if (is.na(occ_file) || is.null(occ_file) || !file.exists(occ_file)) {
      stop("GBIF download unzip succeeded but could not find occurrence.txt/csv inside: ", zip_path)
    }
    
    message("[GBIF] Reading downloaded occurrence file: ", occ_file)
    
    # Read only needed columns (avoid loading huge extra fields)
    needed_cols <- c(
      "gbifID", "occurrenceID",
      "decimalLongitude", "decimalLatitude",
      "eventDate", "year", "countryCode",
      "license",
      "coordinateUncertaintyInMeters",
      "identificationVerificationStatus",
      "issues", "issue",
      "identifiedBy", "dateIdentified",
      # Record type / provenance fields (used for Stage 03+ filtering/diagnostics)
      "basisOfRecord",
      "taxonRank",
      "occurrenceStatus",
      "datasetKey",
      "datasetName",
      "publishingOrgKey",
      "institutionCode",
      "collectionCode"
    )
    
    ct <- readr::cols(
      gbifID = readr::col_character(),
      occurrenceID = readr::col_character(),
      decimalLongitude = readr::col_double(),
      decimalLatitude  = readr::col_double(),
      eventDate = readr::col_character(),
      year = readr::col_integer(),
      countryCode = readr::col_character(),
      license = readr::col_character(),
      coordinateUncertaintyInMeters = readr::col_double(),
      identificationVerificationStatus = readr::col_character(),
      issues = readr::col_character(),
      issue  = readr::col_character(),
      identifiedBy = readr::col_character(),
      dateIdentified = readr::col_character(),
      basisOfRecord = readr::col_character(),
      taxonRank = readr::col_character(),
      occurrenceStatus = readr::col_character(),
      datasetKey = readr::col_character(),
      datasetName = readr::col_character(),
      publishingOrgKey = readr::col_character(),
      institutionCode = readr::col_character(),
      collectionCode = readr::col_character(),
      .default = readr::col_character()
    )
    
    if (grepl("\\.csv$", occ_file, ignore.case = TRUE)) {
      gbif_raw <- readr::read_csv(
        occ_file,
        show_col_types = FALSE,
        progress = TRUE,
        col_select = dplyr::any_of(needed_cols),
        col_types = ct
      )
    } else {
      gbif_raw <- readr::read_tsv(
        occ_file,
        show_col_types = FALSE,
        progress = TRUE,
        col_select = dplyr::any_of(needed_cols),
        col_types = ct
      )
    }
    
    message("GBIF rows read from download file: ", nrow(gbif_raw))
    
    # Harmonise issues column name (download sometimes uses 'issue')
    if (!"issues" %in% names(gbif_raw) && "issue" %in% names(gbif_raw)) {
      gbif_raw$issues <- gbif_raw$issue
    }
    if (!"license" %in% names(gbif_raw)) gbif_raw$license <- NA_character_
    if (!"coordinateUncertaintyInMeters" %in% names(gbif_raw)) gbif_raw$coordinateUncertaintyInMeters <- NA_real_
    if (!"identificationVerificationStatus" %in% names(gbif_raw)) gbif_raw$identificationVerificationStatus <- NA_character_
    if (!"issues" %in% names(gbif_raw)) gbif_raw$issues <- NA_character_
    if (!"identifiedBy" %in% names(gbif_raw)) gbif_raw$identifiedBy <- NA_character_
    if (!"dateIdentified" %in% names(gbif_raw)) gbif_raw$dateIdentified <- NA_character_
    if (!"basisOfRecord" %in% names(gbif_raw)) gbif_raw$basisOfRecord <- NA_character_
    if (!"taxonRank" %in% names(gbif_raw)) gbif_raw$taxonRank <- NA_character_
    if (!"occurrenceStatus" %in% names(gbif_raw)) gbif_raw$occurrenceStatus <- NA_character_
    if (!"datasetKey" %in% names(gbif_raw)) gbif_raw$datasetKey <- NA_character_
    if (!"datasetName" %in% names(gbif_raw)) gbif_raw$datasetName <- NA_character_
    if (!"publishingOrgKey" %in% names(gbif_raw)) gbif_raw$publishingOrgKey <- NA_character_
    if (!"institutionCode" %in% names(gbif_raw)) gbif_raw$institutionCode <- NA_character_
    if (!"collectionCode" %in% names(gbif_raw)) gbif_raw$collectionCode <- NA_character_
    
    # -------------------------------------------------------------------------
    # Basic screening + essential fields (include QA/certainty fields)
    # -------------------------------------------------------------------------
    gbif_clean <- gbif_raw %>%
      transmute(
        source = "GBIF",
        species = species_name,
        gbifID = as.character(gbifID),
        occurrenceID = as.character(occurrenceID),
        lon = decimalLongitude,
        lat = decimalLatitude,
        date = as.character(eventDate),
        year = as.integer(year),
        country = as.character(countryCode),
        
        # Licence fields (no filtering)
        licence_raw = as.character(license),
        licence = lic_normalise_gbif(license),
        licence_expected = !is.na(lic_normalise_gbif(license)) &
          lic_normalise_gbif(license) %in% expected_licences_gbif,
        
        # QA / certainty fields
        coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters),
        identificationVerificationStatus = as.character(identificationVerificationStatus),
        issues = as.character(issues),
        identifiedBy = as.character(identifiedBy),
        dateIdentified = as.character(dateIdentified),
        
        # Record type / provenance fields
        basisOfRecord = as.character(basisOfRecord),
        taxonRank = as.character(taxonRank),
        occurrenceStatus = as.character(occurrenceStatus),
        datasetKey = as.character(datasetKey),
        datasetName = as.character(datasetName),
        publishingOrgKey = as.character(publishingOrgKey),
        institutionCode = as.character(institutionCode),
        collectionCode = as.character(collectionCode)
      ) %>%
      filter(!is.na(lon), !is.na(lat))
    
    message("GBIF rows (download -> screened coords only): ", nrow(gbif_raw), " -> ", nrow(gbif_clean))
    
    # -------------------------------------------------------------------------
    # De-duplication (within GBIF; light-touch)
    # -------------------------------------------------------------------------
    n_before <- nrow(gbif_clean)
    gbif_clean <- gbif_clean %>%
      distinct(gbifID, .keep_all = TRUE) %>%
      distinct(lon, lat, date, .keep_all = TRUE)
    message("GBIF rows (screened -> de-dup): ", n_before, " -> ", nrow(gbif_clean))
    
    # Save output
    write_csv(gbif_clean, gbif_outfile)
    message("Saved GBIF clean file: ", gbif_outfile)
    
    # Mark checkpoint complete
    ckpt$complete <- TRUE
    ckpt$download_status <- "SUCCEEDED"
    ckpt$last_updated <- as.character(Sys.time())
    saveRDS(ckpt, ckpt_file)
    
    attr(gbif_clean, "gbif_status") <- list(state = "complete", method = "download", key = key, expected = total_expected)
    
    # Licence flagging (only if we have records)
    if (nrow(gbif_clean) > 0) {
      unexpected_tbl <- gbif_clean %>%
        mutate(
          licence = as.character(licence),
          licence_raw = as.character(licence_raw),
          licence_expected = !is.na(licence) & licence %in% expected_licences_gbif
        ) %>%
        filter(!licence_expected) %>%
        count(licence_raw, licence, sort = TRUE) %>%
        mutate(
          species = species_name,
          source = "GBIF",
          prop_of_records = n / nrow(gbif_clean)
        )
      
      if (nrow(unexpected_tbl) > 0) {
        message("\n[LICENCE FLAG] GBIF returned licence types outside the expected set for ", species_name, ".")
        message("Expected (normalised): ", paste(expected_licences_gbif, collapse = ", "))
        message("Top unexpected licence entries (see log for full list):")
        print(head(unexpected_tbl, 10), n = 10)
        
        write_unexpected_licence_log(
          species_name = species_name,
          slug = slug,
          source_name = "GBIF",
          unexpected_tbl = unexpected_tbl,
          repo_root = repo_root
        )
      } else {
        message("[LICENCE OK] GBIF: no unexpected licence types detected for ", species_name, ".")
      }
    } else {
      message("[GBIF] No records returned after coordinate screening; skipping licence checks.")
    }
    
    message("GBIF clean: ", nrow(gbif_clean), " records.")
    gbif_clean %>% count(licence, sort = TRUE) %>% print(n = 10)
    
    return(gbif_clean)
  }
  
  # ---------------------------------------------------------------------------
  # SEARCH path (<=100k) — with checkpoint resume, and cap-aware paging
  # ---------------------------------------------------------------------------
  ckpt$mode <- "search"
  ckpt$total_expected <- total_expected
  ckpt$last_updated <- as.character(Sys.time())
  saveRDS(ckpt, ckpt_file)
  
  max_retries <- 5
  retry_base_wait_s <- 10
  
  start <- 0
  all_pages <- list()
  
  # Resume checkpoint if present (older schema used start/all_pages)
  if (!is.null(ckpt$start)) start <- ckpt$start
  if (!is.null(ckpt$all_pages)) all_pages <- ckpt$all_pages
  
  if (start > 0 || length(all_pages) > 0) {
    message("Resuming GBIF search pull from start = ", start,
            " (pages already stored: ", length(all_pages), ").")
  }
  
  repeat {
    Sys.sleep(pause_s)
    
    # Cap-aware paging (prevents start+limit > 100k errors)
    limit_this <- min(page_size, gbif_search_hard_limit - start)
    if (limit_this <= 0) {
      message(
        "\n[GBIF][INCOMPLETE] Hit the occ_search 100k cap for ", species_name, ".\n",
        "This cannot be fixed by waiting/retrying; switch to downloads to retrieve the full dataset.\n"
      )
      ckpt$complete <- FALSE
      ckpt$last_updated <- as.character(Sys.time())
      saveRDS(ckpt, ckpt_file)
      gbif_clean <- empty_gbif_clean()
      attr(gbif_clean, "gbif_status") <- list(state = "capped", method = "search", expected = total_expected)
      return(gbif_clean)
    }
    
    res <- NULL
    for (attempt in seq_len(max_retries)) {
      res <- tryCatch(
        occ_search(
          taxonKey = taxon_key,
          continent = region_scope,
          hasCoordinate = TRUE,
          limit = limit_this,
          start = start
        ),
        error = function(e) e
      )
      
      if (!inherits(res, "error")) break
      
      wait_s <- retry_base_wait_s * attempt
      message("GBIF request failed at start = ", start,
              " (attempt ", attempt, "/", max_retries, "): ",
              conditionMessage(res),
              " | waiting ", wait_s, "s then retrying...")
      Sys.sleep(wait_s)
    }
    
    if (inherits(res, "error")) {
      ckpt$start <- start
      ckpt$all_pages <- all_pages
      ckpt$complete <- FALSE
      ckpt$last_updated <- as.character(Sys.time())
      saveRDS(ckpt, ckpt_file)
      stop(res)
    }
    
    # Total expected from meta if occ_count failed earlier
    if (is.na(total_expected)) total_expected <- res$meta$count
    
    if (length(res$data) == 0) break
    
    all_pages[[length(all_pages) + 1]] <- res$data
    
    ckpt$start <- start
    ckpt$all_pages <- all_pages
    ckpt$total_expected <- total_expected
    ckpt$complete <- FALSE
    ckpt$last_updated <- as.character(Sys.time())
    saveRDS(ckpt, ckpt_file)
    
    pulled_so_far <- start + nrow(res$data)
    message("Pulled ", pulled_so_far, " / ", total_expected, " rows...")
    
    if (pulled_so_far >= total_expected) break
    if (pulled_so_far >= max_records) {
      message("Reached max_records safety limit (", max_records, ").")
      break
    }
    
    start <- start + nrow(res$data)
  }
  
  gbif_raw <- bind_rows(all_pages)
  
  message("GBIF expected rows (query): ", total_expected)
  message("GBIF rows pulled (paged):  ", nrow(gbif_raw))
  
  message("GBIF licence breakdown (RAW pull):")
  gbif_raw %>%
    count(license, sort = TRUE) %>%
    mutate(prop = n / sum(n)) %>%
    print(n = 10)
  
  # -------------------------------------------------------------------------
  # Basic screening + essential fields (include QA/certainty fields)
  # -------------------------------------------------------------------------
  gbif_clean <- gbif_raw %>%
    transmute(
      source = "GBIF",
      species = species_name,
      gbifID = as.character(gbifID),
      occurrenceID = as.character(occurrenceID),
      lon = decimalLongitude,
      lat = decimalLatitude,
      date = as.character(eventDate),
      year = as.integer(year),
      country = as.character(countryCode),
      
      # Licence fields (no filtering)
      licence_raw = as.character(license),
      licence = lic_normalise_gbif(license),
      licence_expected = !is.na(lic_normalise_gbif(license)) &
        lic_normalise_gbif(license) %in% expected_licences_gbif,
      
      # QA / certainty fields
      coordinateUncertaintyInMeters = if ("coordinateUncertaintyInMeters" %in% names(gbif_raw)) {
        as.numeric(coordinateUncertaintyInMeters)
      } else NA_real_,
      
      identificationVerificationStatus = if ("identificationVerificationStatus" %in% names(gbif_raw)) {
        as.character(identificationVerificationStatus)
      } else NA_character_,
      
      issues = if ("issues" %in% names(gbif_raw)) {
        as.character(issues)
      } else NA_character_,
      
      identifiedBy = if ("identifiedBy" %in% names(gbif_raw)) {
        as.character(identifiedBy)
      } else NA_character_,
      
      dateIdentified = if ("dateIdentified" %in% names(gbif_raw)) {
        as.character(dateIdentified)
      } else NA_character_,
      
      # Record type / provenance fields
      basisOfRecord = if ("basisOfRecord" %in% names(gbif_raw)) {
        as.character(basisOfRecord)
      } else NA_character_,
      
      taxonRank = if ("taxonRank" %in% names(gbif_raw)) {
        as.character(taxonRank)
      } else NA_character_,
      
      occurrenceStatus = if ("occurrenceStatus" %in% names(gbif_raw)) {
        as.character(occurrenceStatus)
      } else NA_character_,
      
      datasetKey = if ("datasetKey" %in% names(gbif_raw)) {
        as.character(datasetKey)
      } else NA_character_,
      
      datasetName = if ("datasetName" %in% names(gbif_raw)) {
        as.character(datasetName)
      } else NA_character_,
      
      publishingOrgKey = if ("publishingOrgKey" %in% names(gbif_raw)) {
        as.character(publishingOrgKey)
      } else NA_character_,
      
      institutionCode = if ("institutionCode" %in% names(gbif_raw)) {
        as.character(institutionCode)
      } else NA_character_,
      
      collectionCode = if ("collectionCode" %in% names(gbif_raw)) {
        as.character(collectionCode)
      } else NA_character_
    ) %>%
    filter(!is.na(lon), !is.na(lat))
  
  message("GBIF rows (raw -> screened coords only): ", nrow(gbif_raw), " -> ", nrow(gbif_clean))
  
  # -------------------------------------------------------------------------
  # De-duplication (within GBIF; light-touch)
  # -------------------------------------------------------------------------
  n_before <- nrow(gbif_clean)
  gbif_clean <- gbif_clean %>%
    distinct(gbifID, .keep_all = TRUE) %>%
    distinct(lon, lat, date, .keep_all = TRUE)
  message("GBIF rows (screened -> de-dup): ", n_before, " -> ", nrow(gbif_clean))
  
  # Save output
  write_csv(gbif_clean, gbif_outfile)
  message("Saved GBIF clean file: ", gbif_outfile)
  
  # Mark checkpoint complete only if we retrieved the full dataset (i.e., <=100k)
  ckpt$start <- 0
  ckpt$all_pages <- list()
  ckpt$total_expected <- total_expected
  ckpt$complete <- isTRUE(!is.na(total_expected) && nrow(gbif_raw) >= total_expected && total_expected <= gbif_search_hard_limit)
  ckpt$last_updated <- as.character(Sys.time())
  saveRDS(ckpt, ckpt_file)
  
  # ---------------------------------------------------------------------------
  # Licence flagging
  # ---------------------------------------------------------------------------
  unexpected_tbl <- gbif_clean %>%
    mutate(
      licence = as.character(licence),
      licence_raw = as.character(licence_raw),
      licence_expected = !is.na(licence) & licence %in% expected_licences_gbif
    ) %>%
    filter(!licence_expected) %>%
    count(licence_raw, licence, sort = TRUE) %>%
    mutate(
      species = species_name,
      source = "GBIF",
      prop_of_records = n / nrow(gbif_clean)
    )
  
  if (nrow(unexpected_tbl) > 0) {
    message("\n[LICENCE FLAG] GBIF returned licence types outside the expected set for ", species_name, ".")
    message("Expected (normalised): ", paste(expected_licences_gbif, collapse = ", "))
    message("Top unexpected licence entries (see log for full list):")
    print(head(unexpected_tbl, 10), n = 10)
    
    write_unexpected_licence_log(
      species_name = species_name,
      slug = slug,
      source_name = "GBIF",
      unexpected_tbl = unexpected_tbl,
      repo_root = repo_root
    )
  } else {
    message("[LICENCE OK] GBIF: no unexpected licence types detected for ", species_name, ".")
  }
  
  message("GBIF clean: ", nrow(gbif_clean), " records.")
  gbif_clean %>% count(licence, sort = TRUE) %>% print(n = 10)
  
  attr(gbif_clean, "gbif_status") <- list(
    state = if (isTRUE(ckpt$complete)) "complete" else "incomplete",
    method = "search",
    expected = total_expected
  )
  
  return(gbif_clean)
}

# ==============================================================================
# NBN pull (UK Atlas) ----------------------------------------------------------
# ==============================================================================

pull_nbn_clean <- function(species_name,
                           group_dir = "",
                           species_subdir = FALSE,
                           nbn_email,
                           download_reason_id = 17,
                           expected_licences_nbn = c("OGL", "CC0", "CC-BY", "CC-BY-NC"),
                           use_cache = TRUE,
                           pause_s = 0.25)
{
  
  repo_root <- get_repo_root()
  group_dir <- normalise_group_dir(group_dir)
  nbn_out_root <- if (nzchar(group_dir)) {
    file.path(repo_root, "data", "raw", "nbn", group_dir)
  } else {
    file.path(repo_root, "data", "raw", "nbn")
  }
  
  slug <- slugify_species(species_name)
  
  # Output dirs
  ckpt_root    <- file.path(repo_root, "data", "_checkpoints")
  nbn_ckpt_dir <- file.path(ckpt_root, "nbn")
  
  dir.create(nbn_out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(nbn_ckpt_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Species subfolder (NEW)
  nbn_out_dir <- if (isTRUE(species_subdir)) file.path(nbn_out_root, slug) else nbn_out_root
  dir.create(nbn_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  nbn_outfile   <- file.path(nbn_out_dir,  paste0("nbn_", slug, "_clean.csv"))
  nbn_ckpt_file <- file.path(nbn_ckpt_dir, paste0("nbn_pull_checkpoint_", slug, ".rds"))
  
  # Configure galah to use the UK atlas + provide a download reason
  galah_config(atlas = "United Kingdom", email = nbn_email, verbose = FALSE)
  galah_config(download_reason_id = download_reason_id)
  
  # Normalise NBN licence strings to short codes
  lic_normalise_nbn <- function(x) {
    x_l <- stringr::str_to_lower(x)
    
    dplyr::case_when(
      is.na(x) ~ NA_character_,
      x %in% c("OGL", "CC0", "CC-BY", "CC-BY-NC") ~ x,
      stringr::str_detect(x_l, "open government licence|\\bogl\\b") ~ "OGL",
      stringr::str_detect(x_l, "\\bcc0\\b|publicdomain/zero/1.0") ~ "CC0",
      stringr::str_detect(x_l, "cc[- ]?by[- ]?nc|licenses/by-nc/4.0") ~ "CC-BY-NC",
      stringr::str_detect(x_l, "cc[- ]?by\\b|licenses/by/4.0") ~ "CC-BY",
      TRUE ~ NA_character_
    )
  }
  
  # Columns that MUST exist in cached outputs to be considered "current schema"
  required_cache_cols <- c(
    "coordinateUncertaintyInMeters",
    "identificationVerificationStatus",
    "identifiedBy",
    "coordinatePrecision",
    # Schema alignment / provenance fields (often NA for NBN, but kept for downstream consistency)
    "basisOfRecord",
    "taxonRank",
    "occurrenceStatus",
    "datasetKey",
    "datasetName",
    "publishingOrgKey",
    "institutionCode",
    "collectionCode"
  )
  
  # ---------------------------------------------------------------------------
  # NBN taxon check (kept as-is)
  # ---------------------------------------------------------------------------
  nbn_taxa <- search_taxa(species_name)
  message("NBN taxon search (top hit):")
  nbn_taxa %>%
    select(scientific_name, taxon_concept_id, rank) %>%
    head(1) %>%
    print(n = 1)
  
  # ---------------------------------------------------------------------------
  # NBN occurrence pull (cached if available + schema matches)
  # ---------------------------------------------------------------------------
  use_cached <- isTRUE(use_cache) && file.exists(nbn_outfile)
  
  if (use_cached) {
    nbn_cached <- readr::read_csv(nbn_outfile, show_col_types = FALSE)
    missing_cols <- setdiff(required_cache_cols, names(nbn_cached))
    
    if (length(missing_cols) > 0) {
      message(
        "Found existing NBN clean file but it is missing required columns for the current schema: ",
        paste(missing_cols, collapse = ", "),
        "\nRe-pulling from NBN: ", nbn_outfile
      )
      use_cached <- FALSE
    } else {
      message("Found existing NBN clean file, reading: ", nbn_outfile)
      
      nbn_clean <- nbn_cached %>%
        mutate(
          licence = as.character(licence),
          licence_raw = as.character(licence_raw),
          licence_expected = !is.na(licence) & licence %in% expected_licences_nbn,
          identificationVerificationStatus = as.character(identificationVerificationStatus),
          identifiedBy = as.character(identifiedBy),
          coordinatePrecision = as.character(coordinatePrecision),
          basisOfRecord = as.character(basisOfRecord),
          taxonRank = as.character(taxonRank),
          occurrenceStatus = as.character(occurrenceStatus),
          datasetKey = as.character(datasetKey),
          datasetName = as.character(datasetName),
          publishingOrgKey = as.character(publishingOrgKey),
          institutionCode = as.character(institutionCode),
          collectionCode = as.character(collectionCode)
        )
    }
  }
  
  if (!use_cached) {
    
    max_retries <- 5
    retry_base_wait_s <- 10
    
    # Prefer a checkpoint, but only if it has the new QA fields
    nbn_raw <- NULL
    if (file.exists(nbn_ckpt_file)) {
      message("Found NBN checkpoint, loading: ", nbn_ckpt_file)
      tmp <- readRDS(nbn_ckpt_file)
      
      has_license <- ("dcterms:license" %in% names(tmp)) || ("license" %in% names(tmp))
      has_required <- all(c("recordID", "scientificName", "eventDate", "year",
                            "decimalLatitude", "decimalLongitude") %in% names(tmp))
      has_qa <- all(required_cache_cols %in% names(tmp))
      
      if (has_license && has_required && has_qa) {
        nbn_raw <- tmp
      } else {
        message("Checkpoint exists but is missing QA columns; re-downloading.")
        nbn_raw <- NULL
      }
    }
    
    # If no usable checkpoint, download with retries
    if (is.null(nbn_raw)) {
      
      # Field sets (avoid the known 403 fields: dateIdentified, basisOfRecord, occurrenceStatus)
      nbn_core <- c(
        "recordID",
        "scientificName",
        "eventDate",
        "year",
        "decimalLatitude",
        "decimalLongitude",
        "license"
      )
      
      nbn_qa <- c(
        "coordinateUncertaintyInMeters",
        "coordinatePrecision",
        "identificationVerificationStatus",
        "identifiedBy"
      )
      
      make_select <- function(x) do.call(galah::galah_select, as.list(x))
      
      for (attempt in seq_len(max_retries)) {
        
        Sys.sleep(pause_s)
        
        # 1) Try core + QA in one call
        nbn_raw_try <- tryCatch(
          galah_call() |>
            galah_identify(species_name) |>
            atlas_occurrences(select = make_select(c(nbn_core, nbn_qa))),
          error = function(e) e
        )
        
        if (!inherits(nbn_raw_try, "error")) {
          nbn_raw <- nbn_raw_try
          break
        }
        
        # 2) Fallback: core-only + QA-only (join on recordID)
        message(
          "NBN combined (core+QA) pull failed (attempt ", attempt, "/", max_retries, "): ",
          conditionMessage(nbn_raw_try),
          "\nTrying fallback: core-only + QA-only join..."
        )
        
        core_try <- tryCatch(
          galah_call() |>
            galah_identify(species_name) |>
            atlas_occurrences(select = make_select(nbn_core)),
          error = function(e) e
        )
        
        if (!inherits(core_try, "error")) {
          
          qa_try <- tryCatch(
            galah_call() |>
              galah_identify(species_name) |>
              atlas_occurrences(select = make_select(c("recordID", nbn_qa))),
            error = function(e) e
          )
          
          if (!inherits(qa_try, "error")) {
            qa_try <- qa_try %>% distinct(recordID, .keep_all = TRUE)
            nbn_raw <- core_try %>% left_join(qa_try, by = "recordID")
            break
          } else {
            message("Fallback QA-only pull failed: ", conditionMessage(qa_try))
            # proceed with core only; we'll add QA columns as NA later
            nbn_raw <- core_try
            break
          }
          
        } else {
          wait_s <- retry_base_wait_s * attempt
          message(
            "NBN core pull also failed (attempt ", attempt, "/", max_retries, "): ",
            conditionMessage(core_try),
            " | waiting ", wait_s, "s then retrying..."
          )
          Sys.sleep(wait_s)
        }
      }
      
      if (is.null(nbn_raw)) stop("NBN pull failed after retries for: ", species_name)
      
      # Save raw pull to checkpoint so we don't need to re-download next time
      saveRDS(nbn_raw, nbn_ckpt_file)
      message("Saved NBN checkpoint: ", nbn_ckpt_file)
    }
    
    message("NBN raw rows: ", nrow(nbn_raw))
    
    # -------------------------------------------------------------------------
    # Record check + cleaning (NBN) (UPDATED: include QA/certainty fields)
    # -------------------------------------------------------------------------
    lic_col <- if ("dcterms:license" %in% names(nbn_raw)) "dcterms:license" else "license"
    
    # Ensure QA columns exist even if the fallback ran core-only
    if (!"coordinateUncertaintyInMeters" %in% names(nbn_raw)) nbn_raw$coordinateUncertaintyInMeters <- NA_real_
    if (!"coordinatePrecision" %in% names(nbn_raw))          nbn_raw$coordinatePrecision <- NA_character_
    if (!"identificationVerificationStatus" %in% names(nbn_raw)) nbn_raw$identificationVerificationStatus <- NA_character_
    if (!"identifiedBy" %in% names(nbn_raw))                 nbn_raw$identifiedBy <- NA_character_
    
    # Ensure schema alignment / provenance columns exist (NBN often cannot supply these)
    if (!"basisOfRecord" %in% names(nbn_raw))     nbn_raw$basisOfRecord <- NA_character_
    if (!"taxonRank" %in% names(nbn_raw))         nbn_raw$taxonRank <- NA_character_
    if (!"occurrenceStatus" %in% names(nbn_raw))  nbn_raw$occurrenceStatus <- NA_character_
    if (!"datasetKey" %in% names(nbn_raw))        nbn_raw$datasetKey <- NA_character_
    if (!"datasetName" %in% names(nbn_raw))       nbn_raw$datasetName <- NA_character_
    if (!"publishingOrgKey" %in% names(nbn_raw))  nbn_raw$publishingOrgKey <- NA_character_
    if (!"institutionCode" %in% names(nbn_raw))   nbn_raw$institutionCode <- NA_character_
    if (!"collectionCode" %in% names(nbn_raw))    nbn_raw$collectionCode <- NA_character_
    
    message("NBN licence breakdown (RAW pull):")
    nbn_raw %>%
      count(.data[[lic_col]], sort = TRUE) %>%
      mutate(prop = n / sum(n)) %>%
      print(n = 10)
    
    nbn_clean <- nbn_raw %>%
      transmute(
        source = "NBN",
        species = species_name,
        recordID = as.character(recordID),
        lon = decimalLongitude,
        lat = decimalLatitude,
        date = as.character(eventDate),
        year = as.integer(year),
        
        # Licence fields (no filtering)
        licence_raw = as.character(.data[[lic_col]]),
        licence = lic_normalise_nbn(.data[[lic_col]]),
        licence_expected = !is.na(lic_normalise_nbn(.data[[lic_col]])) &
          lic_normalise_nbn(.data[[lic_col]]) %in% expected_licences_nbn,
        
        # QA / certainty fields (NEW)
        coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters),
        coordinatePrecision = as.character(coordinatePrecision),
        identificationVerificationStatus = as.character(identificationVerificationStatus),
        identifiedBy = as.character(identifiedBy),
        
        # Schema alignment / provenance fields (often NA for NBN)
        basisOfRecord = as.character(basisOfRecord),
        taxonRank = as.character(taxonRank),
        occurrenceStatus = as.character(occurrenceStatus),
        datasetKey = as.character(datasetKey),
        datasetName = as.character(datasetName),
        publishingOrgKey = as.character(publishingOrgKey),
        institutionCode = as.character(institutionCode),
        collectionCode = as.character(collectionCode)
      ) %>%
      filter(!is.na(lon), !is.na(lat))
    
    message("NBN rows (raw -> screened coords only): ", nrow(nbn_raw), " -> ", nrow(nbn_clean))
    
    # -------------------------------------------------------------------------
    # De-duplication (NBN; light-touch)
    # -------------------------------------------------------------------------
    n_before <- nrow(nbn_clean)
    
    nbn_clean <- nbn_clean %>%
      distinct(recordID, .keep_all = TRUE) %>%
      distinct(lon, lat, date, .keep_all = TRUE)
    
    message("NBN rows (screened -> de-dup): ", n_before, " -> ", nrow(nbn_clean))
    
    # Save output (NBN)
    write_csv(nbn_clean, nbn_outfile)
    message("Saved NBN clean file: ", nbn_outfile)
  }
  
  # ---------------------------------------------------------------------------
  # Licence flagging
  # ---------------------------------------------------------------------------
  unexpected_tbl <- nbn_clean %>%
    mutate(
      licence = as.character(licence),
      licence_raw = as.character(licence_raw),
      licence_expected = !is.na(licence) & licence %in% expected_licences_nbn
    ) %>%
    filter(!licence_expected) %>%
    count(licence_raw, licence, sort = TRUE) %>%
    mutate(
      species = species_name,
      source = "NBN",
      prop_of_records = n / nrow(nbn_clean)
    )
  
  if (nrow(unexpected_tbl) > 0) {
    message("\n[LICENCE FLAG] NBN returned licence types outside the expected set for ", species_name, ".")
    message("Expected (normalised): ", paste(expected_licences_nbn, collapse = ", "))
    message("Top unexpected licence entries (see log for full list):")
    print(head(unexpected_tbl, 10), n = 10)
    
    write_unexpected_licence_log(
      species_name = species_name,
      slug = slug,
      source_name = "NBN",
      unexpected_tbl = unexpected_tbl,
      repo_root = repo_root
    )
  } else {
    message("[LICENCE OK] NBN: no unexpected licence types detected for ", species_name, ".")
  }
  
  message("NBN clean: ", nrow(nbn_clean), " records.")
  nbn_clean %>% count(licence, sort = TRUE) %>% print(n = 10)
  
  return(nbn_clean)
}

# ==============================================================================
# Main: pull + save raw outputs for a set of species ----------------------------
# ==============================================================================

pull_raw_occurrences <- function(species_names,
                                 group_dir = "",
                                 species_subdir = FALSE,
                                 region_scope = "EUROPE",
                                 pause_s = 0.25,
                                 page_size = 1000,
                                 max_records = Inf,
                                 nbn_email,
                                 download_reason_id = 17,
                                 expected_licences_gbif = c("CC0_1_0", "CC_BY_4_0", "CC_BY_NC_4_0"),
                                 expected_licences_nbn  = c("OGL", "CC0", "CC-BY", "CC-BY-NC"),
                                 use_cache = TRUE,
                                 gbif_method = c("auto", "search", "download"),
                                 gbif_download_wait = FALSE,
                                 gbif_search_hard_limit = 100000L) {
  
  if (missing(nbn_email) || is.null(nbn_email) || !nzchar(nbn_email)) {
    stop("Please provide nbn_email (the email associated with your NBN Atlas account).")
  }
  
  out <- vector("list", length(species_names))
  names(out) <- species_names
  
  gbif_incomplete <- list()
  
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    message("\n============================================================")
    message("Pulling raw occurrence outputs for: ", sp)
    message("============================================================\n")
    
    gbif_clean <- pull_gbif_clean(
      species_name = sp,
      region_scope = region_scope,
      group_dir = group_dir,
      pause_s = pause_s,
      page_size = page_size,
      max_records = max_records,
      expected_licences_gbif = expected_licences_gbif,
      use_cache = use_cache,
      species_subdir = species_subdir,
      gbif_method = gbif_method,
      gbif_download_wait = gbif_download_wait,
      gbif_search_hard_limit = gbif_search_hard_limit
    )
    
    st <- attr(gbif_clean, "gbif_status")
    if (!is.null(st) && !identical(st$state, "complete")) {
      gbif_incomplete[[sp]] <- st
    }
    
    nbn_clean <- pull_nbn_clean(
      species_name = sp,
      group_dir = group_dir,
      nbn_email = nbn_email,
      download_reason_id = download_reason_id,
      expected_licences_nbn = expected_licences_nbn,
      use_cache = use_cache,
      pause_s = pause_s,
      species_subdir = species_subdir
    )
    
    out[[i]] <- list(gbif_clean = gbif_clean, nbn_clean = nbn_clean)
  }
  
  # ---------------------------------------------------------------------------
  # Final GBIF status warning (e.g., pending downloads)
  # ---------------------------------------------------------------------------
  if (length(gbif_incomplete) > 0) {
    message("\n==================== GBIF WARNING ====================")
    message("Some GBIF pulls are not yet complete.")
    for (nm in names(gbif_incomplete)) {
      s <- gbif_incomplete[[nm]]
      line <- paste0(" - ", nm, ": ", s$state)
      if (!is.null(s$key) && !is.na(s$key) && nzchar(s$key)) line <- paste0(line, " (key=", s$key, ")")
      if (!is.null(s$status) && !is.na(s$status) && nzchar(s$status)) line <- paste0(line, " status=", s$status)
      if (!is.null(s$expected) && !is.na(s$expected)) line <- paste0(line, " expected=", s$expected)
      message(line)
    }
    message("\nRe-run the script later to resume any pending GBIF downloads.")
    message("======================================================\n")
  }
  
  invisible(out)
}
