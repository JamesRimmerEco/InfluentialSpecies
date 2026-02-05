# InfluentialSpecies/R/pull_raw_occurrences.R
#
# InfluentialSpecies — pull raw occurrences (GBIF + NBN Atlas) -------------------
#
# Purpose
#   Pull occurrence data for one or more species from:
#     - GBIF (Europe-wide scope, via rgbif)
#     - NBN Atlas (UK-only, via galah)
#
#   Apply only very basic screening suitable for "raw" outputs:
#     - keep only records with usable coordinates
#     - light-touch, within-source de-duplication only (exact duplicates):
#         GBIF: duplicate gbifID, and exact repeats of (lon, lat, date)
#         NBN : duplicate recordID, and exact repeats of (lon, lat, date)
#
# Outputs
#   For each species we write a per-source "clean" CSV:
#
#   If species_subdir = TRUE:
#     data/raw/gbif/<group_dir>/<slug>/gbif_<slug>_clean.csv
#     data/raw/nbn/<group_dir>/<slug>/nbn_<slug>_clean.csv
#
#   If species_subdir = FALSE:
#     data/raw/gbif/<group_dir>/gbif_<slug>_clean.csv   (or data/raw/gbif/... if group_dir is blank)
#     data/raw/nbn/<group_dir>/nbn_<slug>_clean.csv     (or data/raw/nbn/...  if group_dir is blank)
#
# Checkpoints
#   Checkpoints are stored to support resuming long or asynchronous pulls:
#     data/_checkpoints/gbif/gbif_pull_checkpoint_<slug>.rds
#     data/_checkpoints/nbn/nbn_pull_checkpoint_<slug>.rds
#
#   Optional (recommended on synced/network drives):
#     If you set Sys.setenv(INFLUENTIAL_CHECKPOINT_ROOT = "<local folder>"),
#     checkpoints will be written under that folder instead of inside the repo.
#     This reduces the chance of checkpoint corruption if the repo lives on Google Drive/OneDrive.
#
# How GBIF pulls work (important)
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
#   Reliability improvement (applies to any species, including those <100k):
#   - If occ_search paging repeatedly fails at a particular offset (after retries),
#     the script will (i) retry that page with a smaller page size and, if it still fails,
#     (ii) switch to a GBIF download job for that species so the pipeline can continue.
#
# Typical usage
#   - You usually run this via a wrapper script in /scripts/ (e.g. pull_raw_species_set_*.R).
#   - If no species exceed 100k, one run is enough.
#   - If any species exceed 100k (or if search paging falls back to downloads),
#     you may need to re-run the wrapper script one or more times until the final
#     "GBIF WARNING" list disappears (all downloads completed).
#
# QA / certainty fields included in the "raw clean" outputs
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
# Licence handling
#   We do not filter by licence at this stage.
#   However, we still define "expected" licence sets for each source and, if any other
#   licence types appear, we:
#     (i) flag this clearly to the console, and
#     (ii) write a per-species log file under data/raw/licence_flags/
#
# Notes
#   - If cached *_clean.csv files were created before key fields were added (e.g. QA/provenance),
#     this script treats them as stale and will re-pull.
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
  
  # If credentials are already present in the environment (e.g. set in the console for
  # this session, or provided via .Renviron), do not load a local credentials file that
  # could overwrite them.
  if (nzchar(Sys.getenv("GBIF_USER")) &&
      nzchar(Sys.getenv("GBIF_PWD")) &&
      nzchar(Sys.getenv("GBIF_EMAIL"))) {
    message("GBIF credentials already set in environment; skipping credentials file.")
    return(invisible(TRUE))
  }
  
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
  # A simple, filesystem-safe species "slug" used for filenames and subfolders.
  slug <- str_replace_all(tolower(species_name), "[^a-z0-9]+", "_")
  slug <- str_replace_all(slug, "^_+|_+$", "")
  slug
}

# ---- Helper: make group_dir robust when blank --------------------------------
normalise_group_dir <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) "" else x
}

# ---- Helper: safer checkpoint writing ----------------------------------------
safe_saveRDS <- function(object, file) {
  # Checkpoints are written often; on synced/network drives an interrupted write can leave
  # a corrupt .rds. This writes to a temp file and then replaces the target in one step.
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  
  tmp <- paste0(
    file, ".tmp_",
    format(Sys.time(), "%Y%m%d%H%M%S"),
    "_", sample.int(1e6, 1)
  )
  
  saveRDS(object, tmp)
  
  ok <- file.rename(tmp, file)
  if (!ok) {
    # Fallback path if rename is blocked (e.g. file lock): copy then remove temp.
    ok2 <- file.copy(tmp, file, overwrite = TRUE)
    unlink(tmp)
    if (!ok2) warning("Could not reliably write checkpoint to: ", file)
  }
  
  invisible(TRUE)
}

# ---- Helper: pick a checkpoint root ------------------------------------------
get_checkpoint_root <- function(repo_root) {
  # If INFLUENTIAL_CHECKPOINT_ROOT is set, we write checkpoints there; otherwise use repo/data/_checkpoints.
  ckpt_env <- Sys.getenv("INFLUENTIAL_CHECKPOINT_ROOT")
  if (nzchar(ckpt_env)) return(normalizePath(ckpt_env, winslash = "/", mustWork = FALSE))
  file.path(repo_root, "data", "_checkpoints")
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
                            gbif_download_on_search_error = TRUE,
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
  ckpt_root     <- get_checkpoint_root(repo_root)
  gbif_ckpt_dir <- file.path(ckpt_root, "gbif")
  
  dir.create(gbif_out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(gbif_ckpt_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Species subfolder (optional)
  gbif_out_dir <- if (isTRUE(species_subdir)) file.path(gbif_out_root, slug) else gbif_out_root
  dir.create(gbif_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Clean output file (used as a cache when re-running)
  gbif_outfile <- file.path(gbif_out_dir, paste0("gbif_", slug, "_clean.csv"))
  
  # Checkpoint file (stores download keys + completion state)
  ckpt_file <- file.path(gbif_ckpt_dir, paste0("gbif_pull_checkpoint_", slug, ".rds"))
  
  # Licence normaliser (GBIF uses URL forms in occ_search; download files may already have codes)
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
  
  # Columns that must exist in cached outputs to be considered the current schema
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
  
  # Load / initialise checkpoint (backward compatible with older checkpoint schema)
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
    old <- tryCatch(readRDS(ckpt_file), error = function(e) {
      message("[GBIF] Checkpoint exists but could not be read (will recreate): ", ckpt_file)
      message("       Read error: ", conditionMessage(e))
      NULL
    })
    if (!is.null(old)) ckpt <- utils::modifyList(ckpt, old)
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
      use_cached <- FALSE
    }
  }
  
  # Helper for safe credential check
  have_gbif_creds <- function() {
    nzchar(gbif_user) && nzchar(gbif_pwd) && nzchar(gbif_email)
  }
  
  # ---------------------------------------------------------------------------
  # Decide method: search (<=100k) vs download (>100k) or resume pending download
  # ---------------------------------------------------------------------------
  gbif_method <- match.arg(gbif_method)
  
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
  
  # If a download has already been started for this species, always resume it.
  if (!is.na(ckpt$download_key) && nzchar(ckpt$download_key)) {
    if (!isTRUE(ckpt$complete)) {
      message("[GBIF] Found a pending GBIF download in checkpoint; resuming download mode (key=", ckpt$download_key, ").")
      gbif_method <- "download"
    }
  }
  
  if (gbif_method == "auto") {
    if (!is.na(total_expected) && total_expected > gbif_search_hard_limit) {
      gbif_method <- "download"
    } else {
      gbif_method <- "search"
    }
  }
  
  # If user forced "search" but count suggests it will be capped, switch to downloads (to ensure full data)
  if (gbif_method == "search" && !is.na(total_expected) && total_expected > gbif_search_hard_limit) {
    message("[GBIF] Count exceeds ", gbif_search_hard_limit, ". occ_search() cannot retrieve >100k; switching to downloads.")
    gbif_method <- "download"
  }
  
  # ---------------------------------------------------------------------------
  # DOWNLOAD path (unlimited, async) — submit now, resume on next run
  # ---------------------------------------------------------------------------
  if (gbif_method == "download") {
    
    ckpt$mode <- "download"
    ckpt$total_expected <- total_expected
    ckpt$last_updated <- as.character(Sys.time())
    safe_saveRDS(ckpt, ckpt_file)
    
    if (!have_gbif_creds()) {
      message(
        "\n[GBIF][INCOMPLETE] This species requires GBIF downloads, but GBIF credentials are not available.\n",
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
      
      dl_key <- if (is.list(dl) && "key" %in% names(dl)) dl$key else as.character(dl)
      
      ckpt$download_key <- dl_key
      ckpt$download_status <- "SUBMITTED"
      ckpt$complete <- FALSE
      ckpt$start <- 0
      ckpt$all_pages <- list()
      ckpt$last_updated <- as.character(Sys.time())
      safe_saveRDS(ckpt, ckpt_file)
      
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
    
    key <- ckpt$download_key
    
    if (isTRUE(gbif_download_wait)) {
      message("[GBIF] Waiting for download to finish (key=", key, ") ...")
      rgbif::occ_download_wait(key, status_ping = 30, quiet = FALSE)
    }
    
    meta <- tryCatch(rgbif::occ_download_meta(key), error = function(e) NULL)
    status <- if (!is.null(meta) && !is.null(meta$status)) meta$status else NA_character_
    ckpt$download_status <- status
    ckpt$last_updated <- as.character(Sys.time())
    safe_saveRDS(ckpt, ckpt_file)
    
    if (is.na(status) || !identical(status, "SUCCEEDED")) {
      message(
        "\n[GBIF] Download not ready yet for ", species_name, " (key=", key, ", status=", status, ").\n",
        "Skipping for now — re-run later to resume.\n"
      )
      gbif_clean <- empty_gbif_clean()
      attr(gbif_clean, "gbif_status") <- list(state = "pending_download", method = "download", key = key, status = status, expected = total_expected)
      return(gbif_clean)
    }
    
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
    
    needed_cols <- c(
      "gbifID", "occurrenceID",
      "decimalLongitude", "decimalLatitude",
      "eventDate", "year", "countryCode",
      "license",
      "coordinateUncertaintyInMeters",
      "identificationVerificationStatus",
      "issues", "issue",
      "identifiedBy", "dateIdentified",
      # Record type / provenance fields
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
    
    if (!"issues" %in% names(gbif_raw) && "issue" %in% names(gbif_raw)) gbif_raw$issues <- gbif_raw$issue
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
        
        licence_raw = as.character(license),
        licence = lic_normalise_gbif(license),
        licence_expected = !is.na(lic_normalise_gbif(license)) &
          lic_normalise_gbif(license) %in% expected_licences_gbif,
        
        coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters),
        identificationVerificationStatus = as.character(identificationVerificationStatus),
        issues = as.character(issues),
        identifiedBy = as.character(identifiedBy),
        dateIdentified = as.character(dateIdentified),
        
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
    
    n_before <- nrow(gbif_clean)
    gbif_clean <- gbif_clean %>%
      distinct(gbifID, .keep_all = TRUE) %>%
      distinct(lon, lat, date, .keep_all = TRUE)
    message("GBIF rows (screened -> de-dup): ", n_before, " -> ", nrow(gbif_clean))
    
    write_csv(gbif_clean, gbif_outfile)
    message("Saved GBIF clean file: ", gbif_outfile)
    
    ckpt$complete <- TRUE
    ckpt$download_status <- "SUCCEEDED"
    ckpt$last_updated <- as.character(Sys.time())
    safe_saveRDS(ckpt, ckpt_file)
    
    attr(gbif_clean, "gbif_status") <- list(state = "complete", method = "download", key = key, expected = total_expected)
    
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
    if (nrow(gbif_clean) > 0) gbif_clean %>% count(licence, sort = TRUE) %>% print(n = 10)
    
    return(gbif_clean)
  }
  
  # ---------------------------------------------------------------------------
  # SEARCH path (<=100k) — with checkpoint resume, and cap-aware paging
  # ---------------------------------------------------------------------------
  ckpt$mode <- "search"
  ckpt$total_expected <- total_expected
  ckpt$last_updated <- as.character(Sys.time())
  safe_saveRDS(ckpt, ckpt_file)
  
  max_retries <- 5
  retry_base_wait_s <- 10
  
  start <- 0
  all_pages <- list()
  
  if (!is.null(ckpt$start)) start <- ckpt$start
  if (!is.null(ckpt$all_pages)) all_pages <- ckpt$all_pages
  
  # If the checkpoint claims a non-zero start but has no stored pages, that state is inconsistent.
  # Resetting avoids getting stuck repeatedly resuming from a bad offset without any accumulated data.
  if (isTRUE(start > 0) && length(all_pages) == 0) {
    message("[GBIF] Checkpoint paging state looks inconsistent (start>0 but no stored pages). Resetting paging to start=0.")
    start <- 0
    ckpt$start <- 0
    ckpt$all_pages <- list()
    ckpt$last_updated <- as.character(Sys.time())
    safe_saveRDS(ckpt, ckpt_file)
  }
  
  if (start > 0 || length(all_pages) > 0) {
    message("Resuming GBIF search pull from start = ", start,
            " (pages already stored: ", length(all_pages), ").")
  }
  
  repeat {
    Sys.sleep(pause_s)
    
    limit_this <- min(page_size, gbif_search_hard_limit - start)
    if (limit_this <= 0) {
      message(
        "\n[GBIF][INCOMPLETE] Hit the occ_search 100k cap for ", species_name, ".\n",
        "This cannot be fixed by waiting/retrying; switch to downloads to retrieve the full dataset.\n"
      )
      ckpt$complete <- FALSE
      ckpt$last_updated <- as.character(Sys.time())
      safe_saveRDS(ckpt, ckpt_file)
      gbif_clean <- empty_gbif_clean()
      attr(gbif_clean, "gbif_status") <- list(state = "capped", method = "search", expected = total_expected)
      return(gbif_clean)
    }
    
    res <- NULL
    last_err <- NULL
    
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
      
      last_err <- conditionMessage(res)
      wait_s <- retry_base_wait_s * attempt
      message("GBIF request failed at start = ", start,
              " (attempt ", attempt, "/", max_retries, "): ",
              last_err,
              " | waiting ", wait_s, "s then retrying...")
      Sys.sleep(wait_s)
    }
    
    # If repeated failures occur, retry the same page once using a smaller page size.
    # This keeps the fast page_size for the common case, but increases reliability for problematic pages.
    if (inherits(res, "error")) {
      smaller_limit <- min(300L, limit_this)
      if (smaller_limit < limit_this) {
        message("[GBIF] Retrying the same page with a smaller limit (", smaller_limit, ") to reduce timeout risk...")
        res2 <- tryCatch(
          occ_search(
            taxonKey = taxon_key,
            continent = region_scope,
            hasCoordinate = TRUE,
            limit = smaller_limit,
            start = start
          ),
          error = function(e) e
        )
        if (!inherits(res2, "error")) {
          res <- res2
          limit_this <- smaller_limit
        } else {
          last_err <- conditionMessage(res2)
        }
      }
    }
    
    if (inherits(res, "error")) {
      ckpt$start <- start
      ckpt$all_pages <- all_pages
      ckpt$complete <- FALSE
      ckpt$last_updated <- as.character(Sys.time())
      ckpt$last_error <- last_err
      safe_saveRDS(ckpt, ckpt_file)
      
      # If search paging repeatedly fails, switch to downloads for this species (more robust and resumable).
      if (isTRUE(gbif_download_on_search_error) && have_gbif_creds()) {
        
        message(
          "\n[GBIF][INCOMPLETE] Search paging repeatedly failed for ", species_name, " at start=", start, ".\n",
          "Switching this species to GBIF downloads for reliability (the wrapper can continue).\n",
          "Last error: ", last_err, "\n"
        )
        
        ckpt$mode <- "download"
        ckpt$start <- 0
        ckpt$all_pages <- list()
        ckpt$complete <- FALSE
        ckpt$last_updated <- as.character(Sys.time())
        safe_saveRDS(ckpt, ckpt_file)
        
        # If no download key yet, submit one now.
        if (is.null(ckpt$download_key) || is.na(ckpt$download_key) || !nzchar(ckpt$download_key)) {
          dl <- rgbif::occ_download(
            rgbif::pred_and(
              rgbif::pred("taxonKey", taxon_key),
              rgbif::pred("continent", region_scope),
              rgbif::pred("hasCoordinate", TRUE)
            ),
            user = gbif_user, pwd = gbif_pwd, email = gbif_email
          )
          dl_key <- if (is.list(dl) && "key" %in% names(dl)) dl$key else as.character(dl)
          
          ckpt$download_key <- dl_key
          ckpt$download_status <- "SUBMITTED"
          ckpt$last_updated <- as.character(Sys.time())
          safe_saveRDS(ckpt, ckpt_file)
          
          message("[GBIF] Download submitted (fallback) for ", species_name, " key=", ckpt$download_key)
        } else {
          message("[GBIF] Download already exists in checkpoint; key=", ckpt$download_key)
        }
        
        gbif_clean <- empty_gbif_clean()
        attr(gbif_clean, "gbif_status") <- list(
          state = "pending_download",
          method = "download",
          key = ckpt$download_key,
          expected = total_expected,
          note = "download started after repeated search paging failures"
        )
        return(gbif_clean)
        
      } else {
        message(
          "\n[GBIF][INCOMPLETE] Failed after retries for ", species_name, " at start=", start, ".\n",
          "  Error: ", last_err, "\n",
          "Skipping for now — re-run later to resume.\n"
        )
        
        gbif_clean <- empty_gbif_clean()
        attr(gbif_clean, "gbif_status") <- list(
          state = "error_retry_exhausted",
          method = "search",
          expected = total_expected,
          start = start,
          error = last_err
        )
        return(gbif_clean)
      }
    }
    
    if (is.na(total_expected)) total_expected <- res$meta$count
    if (length(res$data) == 0) break
    
    all_pages[[length(all_pages) + 1]] <- res$data
    
    ckpt$start <- start
    ckpt$all_pages <- all_pages
    ckpt$total_expected <- total_expected
    ckpt$complete <- FALSE
    ckpt$last_updated <- as.character(Sys.time())
    safe_saveRDS(ckpt, ckpt_file)
    
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
  if (nrow(gbif_raw) > 0) {
    gbif_raw %>%
      count(license, sort = TRUE) %>%
      mutate(prop = n / sum(n)) %>%
      print(n = 10)
  } else {
    message("[GBIF] No rows returned in search pull.")
  }
  
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
      
      licence_raw = as.character(license),
      licence = lic_normalise_gbif(license),
      licence_expected = !is.na(lic_normalise_gbif(license)) &
        lic_normalise_gbif(license) %in% expected_licences_gbif,
      
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
  
  n_before <- nrow(gbif_clean)
  gbif_clean <- gbif_clean %>%
    distinct(gbifID, .keep_all = TRUE) %>%
    distinct(lon, lat, date, .keep_all = TRUE)
  message("GBIF rows (screened -> de-dup): ", n_before, " -> ", nrow(gbif_clean))
  
  write_csv(gbif_clean, gbif_outfile)
  message("Saved GBIF clean file: ", gbif_outfile)
  
  ckpt$start <- 0
  ckpt$all_pages <- list()
  ckpt$total_expected <- total_expected
  ckpt$complete <- isTRUE(!is.na(total_expected) && nrow(gbif_raw) >= total_expected && total_expected <= gbif_search_hard_limit)
  ckpt$last_updated <- as.character(Sys.time())
  safe_saveRDS(ckpt, ckpt_file)
  
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
  if (nrow(gbif_clean) > 0) gbif_clean %>% count(licence, sort = TRUE) %>% print(n = 10)
  
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
  ckpt_root    <- get_checkpoint_root(repo_root)
  nbn_ckpt_dir <- file.path(ckpt_root, "nbn")
  
  dir.create(nbn_out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(nbn_ckpt_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Species subfolder (optional)
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
  
  # Columns that must exist in cached outputs to be considered the current schema
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
  
  # Always return a tibble with the expected schema (even if 0 rows)
  empty_nbn_clean <- function() {
    tibble::tibble(
      source = character(),
      species = character(),
      recordID = character(),
      lon = numeric(),
      lat = numeric(),
      date = character(),
      year = integer(),
      licence_raw = character(),
      licence = character(),
      licence_expected = logical(),
      coordinateUncertaintyInMeters = numeric(),
      coordinatePrecision = character(),
      identificationVerificationStatus = character(),
      identifiedBy = character(),
      # Schema alignment / provenance fields (often NA for NBN)
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
      
      if (nrow(nbn_clean) == 0) {
        message("[NBN] Cached file is empty; skipping licence checks.")
        message("NBN clean: 0 records.")
        return(nbn_clean)
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # NBN taxon guard
  #   We proceed only if NBN taxonomy contains an exact, species-rank match.
  #   This reduces the chance of pulling near-matches or surprising synonyms.
  #   If there is no exact match, we write an empty (schema-correct) CSV and continue.
  # ---------------------------------------------------------------------------
  if (!use_cached) {
    
    nbn_taxa <- search_taxa(species_name)
    message("NBN taxon search (top hit):")
    
    if (inherits(nbn_taxa, "data.frame")) {
      nbn_taxa %>%
        dplyr::select(dplyr::any_of(c("scientific_name", "scientificName",
                                      "taxon_concept_id", "taxonConceptId",
                                      "rank"))) %>%
        head(1) %>%
        print(n = 1)
    } else {
      print(utils::head(nbn_taxa, 1))
    }
    
    nbn_taxa2 <- nbn_taxa
    
    if (inherits(nbn_taxa2, "data.frame")) {
      
      if (!"scientific_name" %in% names(nbn_taxa2)) {
        if ("scientificName" %in% names(nbn_taxa2)) {
          nbn_taxa2$scientific_name <- nbn_taxa2$scientificName
        } else {
          nbn_taxa2$scientific_name <- NA_character_
        }
      }
      
      if (!"rank" %in% names(nbn_taxa2)) {
        nbn_taxa2$rank <- NA_character_
      }
      
      nbn_exact <- nbn_taxa2 %>%
        filter(
          !is.na(scientific_name),
          tolower(scientific_name) == tolower(species_name),
          !is.na(rank),
          tolower(rank) == "species"
        )
      
      if (nrow(nbn_exact) == 0) {
        message(
          "[NBN] No exact species match for '", species_name, "' in NBN taxonomy.\n",
          "      (Non-UK taxon, synonym/spelling difference, or absent from NBN.)\n",
          "      Skipping NBN pull and writing an empty output so the pipeline can continue."
        )
        
        nbn_clean <- empty_nbn_clean()
        readr::write_csv(nbn_clean, nbn_outfile)
        message("Saved NBN clean file (EMPTY): ", nbn_outfile)
        return(nbn_clean)
      }
    } else {
      message("[NBN] Taxon table format unexpected; proceeding to attempt pull.")
    }
  }
  
  max_retries <- 5
  retry_base_wait_s <- 10
  
  # Prefer a checkpoint, but only if it has the required fields for the current schema
  nbn_raw <- NULL
  if (file.exists(nbn_ckpt_file)) {
    message("Found NBN checkpoint, loading: ", nbn_ckpt_file)
    tmp <- tryCatch(readRDS(nbn_ckpt_file), error = function(e) {
      message("[NBN] Checkpoint exists but could not be read (will re-download): ", nbn_ckpt_file)
      message("      Read error: ", conditionMessage(e))
      NULL
    })
    
    if (!is.null(tmp)) {
      has_license <- ("dcterms:license" %in% names(tmp)) || ("license" %in% names(tmp))
      has_required <- all(c("recordID", "scientificName", "eventDate", "year",
                            "decimalLatitude", "decimalLongitude") %in% names(tmp))
      has_qa <- all(required_cache_cols %in% names(tmp))
      
      if (has_license && has_required && has_qa) {
        nbn_raw <- tmp
      } else {
        message("Checkpoint exists but is missing required columns; re-downloading.")
        nbn_raw <- NULL
      }
    }
  }
  
  # If no usable checkpoint, download with retries
  if (is.null(nbn_raw)) {
    
    # Field sets (avoid known 403 fields: dateIdentified, basisOfRecord, occurrenceStatus)
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
    
    safe_saveRDS(nbn_raw, nbn_ckpt_file)
    message("Saved NBN checkpoint: ", nbn_ckpt_file)
  }
  
  message("NBN raw rows: ", nrow(nbn_raw))
  
  lic_col <- if ("dcterms:license" %in% names(nbn_raw)) "dcterms:license" else "license"
  
  if (!"coordinateUncertaintyInMeters" %in% names(nbn_raw)) nbn_raw$coordinateUncertaintyInMeters <- NA_real_
  if (!"coordinatePrecision" %in% names(nbn_raw))          nbn_raw$coordinatePrecision <- NA_character_
  if (!"identificationVerificationStatus" %in% names(nbn_raw)) nbn_raw$identificationVerificationStatus <- NA_character_
  if (!"identifiedBy" %in% names(nbn_raw))                 nbn_raw$identifiedBy <- NA_character_
  
  if (!"basisOfRecord" %in% names(nbn_raw))     nbn_raw$basisOfRecord <- NA_character_
  if (!"taxonRank" %in% names(nbn_raw))         nbn_raw$taxonRank <- NA_character_
  if (!"occurrenceStatus" %in% names(nbn_raw))  nbn_raw$occurrenceStatus <- NA_character_
  if (!"datasetKey" %in% names(nbn_raw))        nbn_raw$datasetKey <- NA_character_
  if (!"datasetName" %in% names(nbn_raw))       nbn_raw$datasetName <- NA_character_
  if (!"publishingOrgKey" %in% names(nbn_raw))  nbn_raw$publishingOrgKey <- NA_character_
  if (!"institutionCode" %in% names(nbn_raw))   nbn_raw$institutionCode <- NA_character_
  if (!"collectionCode" %in% names(nbn_raw))    nbn_raw$collectionCode <- NA_character_
  
  if (nrow(nbn_raw) > 0) {
    message("NBN licence breakdown (RAW pull):")
    nbn_raw %>%
      count(.data[[lic_col]], sort = TRUE) %>%
      mutate(prop = n / sum(n)) %>%
      print(n = 10)
  } else {
    message("[NBN] No rows returned (0 UK records is plausible for non-native taxa).")
  }
  
  nbn_clean <- nbn_raw %>%
    transmute(
      source = "NBN",
      species = species_name,
      recordID = as.character(recordID),
      lon = decimalLongitude,
      lat = decimalLatitude,
      date = as.character(eventDate),
      year = as.integer(year),
      
      licence_raw = as.character(.data[[lic_col]]),
      licence = lic_normalise_nbn(.data[[lic_col]]),
      licence_expected = !is.na(lic_normalise_nbn(.data[[lic_col]])) &
        lic_normalise_nbn(.data[[lic_col]]) %in% expected_licences_nbn,
      
      coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters),
      coordinatePrecision = as.character(coordinatePrecision),
      identificationVerificationStatus = as.character(identificationVerificationStatus),
      identifiedBy = as.character(identifiedBy),
      
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
  
  n_before <- nrow(nbn_clean)
  
  nbn_clean <- nbn_clean %>%
    distinct(recordID, .keep_all = TRUE) %>%
    distinct(lon, lat, date, .keep_all = TRUE)
  
  message("NBN rows (screened -> de-dup): ", n_before, " -> ", nrow(nbn_clean))
  
  write_csv(nbn_clean, nbn_outfile)
  message("Saved NBN clean file: ", nbn_outfile)
  
  if (nrow(nbn_clean) == 0) {
    message("[NBN] No records after coordinate screening; skipping licence checks.")
    message("NBN clean: 0 records.")
    return(nbn_clean)
  }
  
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
                                 gbif_search_hard_limit = 100000L,
                                 gbif_download_on_search_error = TRUE) {
  
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
      gbif_search_hard_limit = gbif_search_hard_limit,
      gbif_download_on_search_error = gbif_download_on_search_error
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
  # Final GBIF status warning (e.g., pending downloads, retry-exhausted search)
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
      if (!is.null(s$start) && !is.na(s$start)) line <- paste0(line, " start=", s$start)
      if (!is.null(s$note) && !is.na(s$note) && nzchar(s$note)) line <- paste0(line, " (", s$note, ")")
      message(line)
    }
    message("\nRe-run the script later to resume any pending GBIF downloads.")
    message("======================================================\n")
  }
  
  invisible(out)
}
