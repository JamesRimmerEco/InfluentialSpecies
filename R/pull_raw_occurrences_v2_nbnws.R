# InfluentialSpecies/R/pull_raw_occurrences_v2_nbnws.R
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
#     <checkpoint_root>/gbif/gbif_pull_checkpoint_<slug>.rds
#     <checkpoint_root>/nbn/nbn_state_<slug>.rds
#
#   Optional (recommended on synced/network drives):
#     If you set Sys.setenv(INFLUENTIAL_CHECKPOINT_ROOT = "<local folder>"),
#     checkpoints will be written under that folder instead of inside the repo.
#     This reduces the chance of checkpoint corruption if the repo lives on Google Drive/OneDrive.
#
# GBIF work files (disk safety)
#   GBIF downloads can be huge (multi-GB zips + much larger extracted files). By default we:
#     - download zips into INFLUENTIAL_GBIF_WORK_ROOT/gbif_zips
#     - extract into INFLUENTIAL_GBIF_WORK_ROOT/gbif_unzip
#     - delete both the zip and extracted folder once the clean CSV is written successfully
#
#   If INFLUENTIAL_GBIF_WORK_ROOT is not set, we fall back to the checkpoint root.
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
# Licence handling
#   We do not filter by licence at this stage.
#   However, we still define "expected" licence sets for each source and, if any other
#   licence types appear, we:
#     (i) flag this clearly to the console, and
#     (ii) write a per-species log file under data/raw/licence_flags/
# ------------------------------------------------------------------------------

# Guard: a broken na.print option can crash printing (error: invalid 'na.print' specification)
opt_na_print <- getOption("na.print")
if (!is.character(opt_na_print) || length(opt_na_print) != 1L || is.na(opt_na_print)) {
  options(na.print = "NA")
}

# Helper: ensure na.print is safe before any console printing of tibbles/data.frames.
# This is defensive against sessions where options(na.print=...) has been set to an invalid value.
ensure_safe_na_print <- function() {
  opt <- getOption("na.print")
  if (!is.character(opt) || length(opt) != 1L || is.na(opt)) options(na.print = "NA")
  invisible(TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(lubridate)
  library(tibble)
  library(rgbif)
  library(galah)
  library(jsonlite)
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

# ---- Helper: GBIF work folders (zips + extraction) ---------------------------
get_gbif_work_root <- function(repo_root) {
  # This is where big GBIF artefacts live briefly (zip + extraction), and are cleaned after success.
  root <- Sys.getenv("INFLUENTIAL_GBIF_WORK_ROOT")
  if (nzchar(root)) return(normalizePath(root, winslash = "/", mustWork = FALSE))
  get_checkpoint_root(repo_root)
}

gbif_work_dirs <- function(repo_root) {
  root <- get_gbif_work_root(repo_root)
  zip_dir <- file.path(root, "gbif_zips")
  unzip_root <- file.path(root, "gbif_unzip")
  dir.create(zip_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(unzip_root, recursive = TRUE, showWarnings = FALSE)
  list(root = root, zip_dir = zip_dir, unzip_root = unzip_root)
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

# ---- Helper: robust read of GBIF download occurrence file --------------------
read_gbif_download_occurrence <- function(occ_file, needed_cols) {
  
  # GBIF downloads are Darwin Core; occurrence.txt is normally tab-separated with many columns.
  # Reading the full file into memory is expensive; we pull only the columns we actually use.
  delim <- if (grepl("\\.csv$", occ_file, ignore.case = TRUE)) "," else "\t"
  
  hdr <- readLines(occ_file, n = 1, warn = FALSE, encoding = "UTF-8")
  if (length(hdr) == 0) stop("GBIF occurrence file appears to be empty: ", occ_file)
  
  hdr_cols <- strsplit(hdr, split = delim, fixed = TRUE)[[1]]
  idx <- match(needed_cols, hdr_cols)
  keep <- which(!is.na(idx))
  if (length(keep) == 0) stop("Could not find any expected columns in GBIF download: ", occ_file)
  
  # Prefer data.table::fread if available (it is much more robust for very large tab files).
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(
      occ_file,
      sep = delim,
      select = idx[keep],
      showProgress = TRUE,
      quote = "",
      encoding = "UTF-8"
    )
    df <- as.data.frame(dt)
  } else {
    # Fallback to readr (works, but more fragile for very large files)
    df <- readr::read_delim(
      occ_file,
      delim = delim,
      show_col_types = FALSE,
      progress = TRUE,
      col_select = dplyr::all_of(hdr_cols[idx[keep]]),
      name_repair = "minimal"
    )
    df <- as.data.frame(df)
  }
  
  # Standardise types we care about
  if ("decimalLongitude" %in% names(df)) df$decimalLongitude <- as.numeric(df$decimalLongitude)
  if ("decimalLatitude"  %in% names(df)) df$decimalLatitude  <- as.numeric(df$decimalLatitude)
  if ("year" %in% names(df)) df$year <- suppressWarnings(as.integer(df$year))
  if ("coordinateUncertaintyInMeters" %in% names(df)) df$coordinateUncertaintyInMeters <- as.numeric(df$coordinateUncertaintyInMeters)
  
  df
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
                            cleanup_gbif_work_files = TRUE,
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
  clean_species_name <- function(x) {
    # Defensive normalisation: trim, collapse whitespace, and convert non-breaking space.
    x <- as.character(x)
    x <- gsub("\u00A0", " ", x, fixed = TRUE)
    x <- trimws(x)
    x <- gsub("\\s+", " ", x)
    x
  }
  
  sp_clean <- clean_species_name(species_name)
  bb <- tryCatch(
    rgbif::name_backbone(name = sp_clean, kingdom = "Animalia", rank = "species"),
    error = function(e) e
  )
  
  if (inherits(bb, "error")) {
    msg <- conditionMessage(bb)
    message("GBIF match:  (usageKey=, matchType=)")
    message(
      "\n[GBIF][INCOMPLETE] Could not resolve GBIF taxonKey for ", species_name, ".\n",
      "  Error: ", msg, "\n",
      "Skipping GBIF pull and returning an empty output so the pipeline can continue.\n"
    )
    gbif_clean <- empty_gbif_clean()
    attr(gbif_clean, "gbif_status") <- list(state = "taxon_unresolved", method = "none", expected = NA_integer_, error = msg)
    return(gbif_clean)
  }
  
  taxon_key <- if ("usageKey" %in% names(bb)) bb$usageKey[1] else NA
  message("GBIF match: ", bb$scientificName, " (usageKey=", taxon_key, ", matchType=", bb$matchType, ")")
  
  if (is.null(taxon_key) || length(taxon_key) != 1L || is.na(taxon_key) || !nzchar(as.character(taxon_key))) {
    note <- if ("note" %in% names(bb)) as.character(bb$note[1]) else NA_character_
    message(
      "\n[GBIF][INCOMPLETE] Could not resolve a unique GBIF taxonKey for ", species_name, ".\n",
      "Skipping GBIF pull and returning an empty output so the pipeline can continue.\n"
    )
    if (!is.na(note) && nzchar(note)) message("  note: ", note)
    gbif_clean <- empty_gbif_clean()
    attr(gbif_clean, "gbif_status") <- list(
      state = "taxon_unresolved",
      method = "none",
      expected = NA_integer_,
      matchType = as.character(bb$matchType[1]),
      note = note
    )
    return(gbif_clean)
  }
  
  taxon_key <- as.integer(taxon_key)
  
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
      
      dl <- tryCatch(
        rgbif::occ_download(
          rgbif::pred_and(
            rgbif::pred("taxonKey", taxon_key),
            rgbif::pred("continent", region_scope),
            rgbif::pred("hasCoordinate", TRUE)
          ),
          user = gbif_user, pwd = gbif_pwd, email = gbif_email
        ),
        error = function(e) e
      )
      
      if (inherits(dl, "error")) {
        msg <- conditionMessage(dl)
        gbif_clean <- empty_gbif_clean()
        attr(gbif_clean, "gbif_status") <- list(state = "download_submit_failed", method = "download", expected = total_expected, error = msg)
        stop(msg)
      }
      
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
    
    # Work dirs for zip + extraction (cleaned after success)
    wd <- gbif_work_dirs(repo_root)
    
    zip_path <- rgbif::occ_download_get(key, path = wd$zip_dir, overwrite = TRUE)
    if (is.list(zip_path) && "path" %in% names(zip_path)) zip_path <- zip_path$path
    
    unzip_dir <- file.path(wd$unzip_root, paste0("gbif_dwc_", key))
    dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Always clear extraction folders; zip is removed only after a clean success.
    on.exit({
      if (isTRUE(cleanup_gbif_work_files)) {
        unlink(unzip_dir, recursive = TRUE, force = TRUE)
      }
    }, add = TRUE)
    
    utils::unzip(zip_path, exdir = unzip_dir)
    
    occ_file <- list.files(
      unzip_dir,
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
    
    gbif_raw <- read_gbif_download_occurrence(occ_file, needed_cols)
    
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
        lon = as.numeric(decimalLongitude),
        lat = as.numeric(decimalLatitude),
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
    
    # Clean up big GBIF artefacts as soon as we have the clean CSV on disk.
    if (isTRUE(cleanup_gbif_work_files)) {
      unlink(zip_path, force = TRUE)
    }
    
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
        
        tryCatch(
          {
            ensure_safe_na_print()
            print(dplyr::slice_head(unexpected_tbl, n = 10), n = 10)
          },
          error = function(e) {
            ensure_safe_na_print()
            message("[GBIF] NOTE: could not print unexpected licence table (", conditionMessage(e), "). Continuing.")
          }
        )
        
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
    if (nrow(gbif_clean) > 0) {
      tryCatch(
        {
          ensure_safe_na_print()
          gbif_clean %>% count(licence, sort = TRUE) %>% print(n = 10)
        },
        error = function(e) {
          ensure_safe_na_print()
          message("[GBIF] NOTE: could not print final licence table (", conditionMessage(e), "). Continuing.")
        }
      )
    }
    
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
    tryCatch(
      {
        ensure_safe_na_print()
        gbif_raw %>%
          count(license, sort = TRUE) %>%
          mutate(prop = n / sum(n)) %>%
          print(n = 10)
      },
      error = function(e) {
        ensure_safe_na_print()
        message("[GBIF] NOTE: could not print licence breakdown (", conditionMessage(e), "). Continuing.")
      }
    )
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
      
      tryCatch(
        {
          ensure_safe_na_print()
          print(dplyr::slice_head(unexpected_tbl, n = 10), n = 10)
        },
        error = function(e) {
          ensure_safe_na_print()
          message("[GBIF] NOTE: could not print unexpected licence table (", conditionMessage(e), "). Continuing.")
        }
      )
      
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
  if (nrow(gbif_clean) > 0) {
    tryCatch(
      {
        ensure_safe_na_print()
        gbif_clean %>% count(licence, sort = TRUE) %>% print(n = 10)
      },
      error = function(e) {
        ensure_safe_na_print()
        message("[GBIF] NOTE: could not print final licence table (", conditionMessage(e), "). Continuing.")
      }
    )
  }
  
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
                           pause_s = 0.25,
                           nbn_download_timeout_s = 3600L) {
  
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
  
  nbn_outfile <- file.path(nbn_out_dir,  paste0("nbn_", slug, "_clean.csv"))
  
  # NBN completion state (small, avoids the "empty CSV looks done forever" problem)
  nbn_state_file <- file.path(nbn_ckpt_dir, paste0("nbn_state_", slug, ".rds"))
  nbn_state <- list(complete = FALSE, last_updated = NA_character_, note = NA_character_, last_error = NA_character_,
                    totalRecords = NA_integer_, guid = NA_character_)
  if (file.exists(nbn_state_file)) {
    tmp <- tryCatch(readRDS(nbn_state_file), error = function(e) NULL)
    if (!is.null(tmp) && is.list(tmp)) nbn_state <- utils::modifyList(nbn_state, tmp)
  }
  
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
  
  # ---- NBN web services fallback (species-ws + records-ws) -------------------
  # galah relies on the NBN "species-ws" taxonomy service for name resolution.
  # For some taxa, that service returns non-rectangular fields (e.g. synonymComplete
  # as a list), which can trigger a deterministic parsing error inside galah.
  #
  # When this happens we:
  #   1) resolve the taxon GUID ourselves via species-ws; then
  #   2) pull occurrences via records-ws, bypassing galah’s taxonomy parser.
  #
  # This keeps output folders, filenames, and schemas unchanged.
  
  nbn_species_ws_search <- function(sp) {
    u <- paste0(
      "https://species-ws.nbnatlas.org/search?q=",
      utils::URLencode(sp, reserved = TRUE),
      "&fq=idxtype:TAXON&pageSize=50"
    )
    
    raw <- tryCatch(
      jsonlite::fromJSON(u, simplifyVector = FALSE, simplifyDataFrame = FALSE),
      error = function(e) e
    )
    if (inherits(raw, "error")) stop("NBN species-ws search failed: ", conditionMessage(raw))
    
    res <- raw$searchResults$results
    if (is.null(res) || length(res) == 0) return(data.frame())
    
    get1 <- function(x, nm, alt = NULL) {
      v <- x[[nm]]
      if ((is.null(v) || length(v) == 0) && !is.null(alt)) v <- x[[alt]]
      if (is.null(v) || length(v) == 0) return(NA_character_)
      as.character(v[[1]])
    }
    
    data.frame(
      scientificName   = vapply(res, get1, character(1), nm = "scientificName", alt = "name"),
      rank             = vapply(res, get1, character(1), nm = "rank"),
      taxonomicStatus  = vapply(res, get1, character(1), nm = "taxonomicStatus"),
      guid             = vapply(res, get1, character(1), nm = "guid"),
      occurrenceCount  = suppressWarnings(as.integer(vapply(res, get1, character(1), nm = "occurrenceCount"))),
      stringsAsFactors = FALSE
    )
  }
  
  nbn_pick_guid <- function(res, sp) {
    if (!is.data.frame(res) || nrow(res) == 0) return(NA_character_)
    
    # Standardise a few key columns we expect to exist (but keep this defensive).
    if (!"scientificName" %in% names(res) && "name" %in% names(res)) res$scientificName <- res$name
    if (!"rank" %in% names(res)) res$rank <- NA_character_
    if (!"taxonomicStatus" %in% names(res)) res$taxonomicStatus <- NA_character_
    if (!"guid" %in% names(res)) res$guid <- NA_character_
    if (!"occurrenceCount" %in% names(res)) res$occurrenceCount <- NA_integer_
    
    res2 <- res %>%
      mutate(
        scientificName2 = tolower(as.character(scientificName)),
        rank2 = tolower(as.character(rank)),
        status2 = tolower(as.character(taxonomicStatus)),
        occ_n = suppressWarnings(as.integer(occurrenceCount))
      )
    
    # First choice: exact binomial, species-rank, accepted.
    hit <- res2 %>%
      filter(
        !is.na(scientificName2),
        scientificName2 == tolower(sp),
        !is.na(rank2),
        rank2 == "species"
      ) %>%
      mutate(is_accepted = !is.na(status2) & status2 == "accepted") %>%
      arrange(desc(is_accepted), desc(occ_n)) %>%
      slice_head(n = 1)
    
    if (nrow(hit) == 0) return(NA_character_)
    guid <- as.character(hit$guid[1])
    if (!nzchar(guid)) return(NA_character_)
    guid
  }
  
  nbn_build_query <- function(params) {
    paste(
      paste0(
        names(params), "=",
        vapply(params, function(x) utils::URLencode(as.character(x), reserved = TRUE), character(1))
      ),
      collapse = "&"
    )
  }
  
  nbn_is_zip_file <- function(path) {
    if (!file.exists(path)) return(FALSE)
    sig <- tryCatch(readBin(path, "raw", n = 2), error = function(e) raw(0))
    length(sig) == 2 && identical(sig, as.raw(c(0x50, 0x4B))) # "PK"
  }
  
  nbn_read_download_file <- function(path) {
    if (!file.exists(path)) stop("NBN download file not found: ", path)
    
    if (requireNamespace("data.table", quietly = TRUE)) {
      df <- data.table::fread(path, encoding = "UTF-8", showProgress = TRUE)
      as.data.frame(df)
    } else {
      readr::read_delim(
        path,
        delim = if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) "\t" else ",",
        show_col_types = FALSE,
        progress = TRUE,
        name_repair = "minimal"
      ) %>%
        as.data.frame()
    }
  }
  
  # Choose the occurrence export file from an unzipped download.
  # NBN downloads can include small helper CSVs (e.g. headings/citation); we need the actual occurrence table.
  nbn_read_header_line <- function(path) {
    tryCatch(readLines(path, n = 1, warn = FALSE, encoding = "UTF-8"), error = function(e) "")
  }
  
  nbn_header_has_coords <- function(hdr) {
    if (!nzchar(hdr)) return(FALSE)
    h <- tolower(hdr)
    grepl("decimallatitude|decimallongitude|gridreference|grid_ref|easting|northing", h)
  }
  
  nbn_choose_occurrence_data_file <- function(files) {
    files <- files[file.exists(files)]
    if (length(files) == 0) return(NA_character_)
    
    base <- tolower(basename(files))
    drop <- base %in% c("headings.csv", "heading.csv", "citation.csv", "citations.csv", "readme.txt", "readme.csv", "metadata.csv")
    keep_files <- files[!drop]
    if (length(keep_files) == 0) keep_files <- files
    
    hdrs <- vapply(keep_files, nbn_read_header_line, character(1))
    has_coords <- vapply(hdrs, nbn_header_has_coords, logical(1))
    sizes <- suppressWarnings(as.numeric(file.info(keep_files)$size))
    sizes[is.na(sizes)] <- 0
    
    if (any(has_coords)) {
      cand <- keep_files[has_coords]
      cand_sizes <- sizes[has_coords]
      return(cand[which.max(cand_sizes)][1])
    }
    
    keep_files[which.max(sizes)][1]
  }
  
  nbn_standardise_ws_raw <- function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) {
      out <- data.frame(
        recordID = character(),
        scientificName = character(),
        eventDate = character(),
        year = integer(),
        decimalLatitude = numeric(),
        decimalLongitude = numeric(),
        license = character(),
        coordinateUncertaintyInMeters = numeric(),
        coordinatePrecision = character(),
        identificationVerificationStatus = character(),
        identifiedBy = character(),
        basisOfRecord = character(),
        taxonRank = character(),
        occurrenceStatus = character(),
        datasetKey = character(),
        datasetName = character(),
        publishingOrgKey = character(),
        institutionCode = character(),
        collectionCode = character(),
        stringsAsFactors = FALSE
      )
      return(out)
    }
    
    # Case-insensitive column picker (keeps the original column name)
    pick_col <- function(candidates) {
      nms <- names(df)
      map <- stats::setNames(nms, tolower(nms))
      cand_l <- tolower(candidates)
      hit_l <- cand_l[cand_l %in% names(map)]
      if (length(hit_l) == 0) return(NULL)
      unname(map[hit_l[1]])
    }
    
    col_record <- pick_col(c("recordID", "recordId", "record_uuid", "uuid", "id"))
    col_sci    <- pick_col(c("scientificName", "scientific_name", "taxon_name", "species"))
    col_date   <- pick_col(c("eventDate", "event_date", "eventdate", "date", "occurrence_date"))
    col_year   <- pick_col(c("year", "eventYear"))
    col_lat    <- pick_col(c("decimalLatitude", "decimal_latitude", "latitude", "lat"))
    col_lon    <- pick_col(c("decimalLongitude", "decimal_longitude", "longitude", "lon", "lng"))
    col_lic    <- pick_col(c("license", "licence", "dcterms:license", "dcterms.license"))
    
    col_cuim   <- pick_col(c("coordinateUncertaintyInMeters", "coordinate_uncertainty_in_meters", "coord_uncertainty_m"))
    col_cp     <- pick_col(c("coordinatePrecision", "coordinate_precision"))
    col_iv     <- pick_col(c("identificationVerificationStatus", "identification_verification_status", "verificationstatus", "verification_status"))
    col_idby   <- pick_col(c("identifiedBy", "identified_by"))
    
    # Some records-ws exports include a point field like "lat,lon" (e.g. point-1km / point-100m).
    col_point  <- pick_col(c("point00001", "point0001", "point001", "point01", "point1", "point", "point_1km", "point_100m", "point_10m"))
    col_grid   <- pick_col(c("gridReference", "grid_reference", "grid_ref"))
    
    # Provenance-ish fields (often absent from NBN downloads; we keep them for schema alignment)
    col_basis  <- pick_col(c("basisOfRecord", "basis_of_record"))
    col_rank   <- pick_col(c("taxonRank", "taxon_rank", "rank"))
    col_occst  <- pick_col(c("occurrenceStatus", "occurrence_status"))
    col_dk     <- pick_col(c("datasetKey", "dataset_key"))
    col_dn     <- pick_col(c("datasetName", "dataset_name"))
    col_pok    <- pick_col(c("publishingOrgKey", "publishing_org_key"))
    col_inst   <- pick_col(c("institutionCode", "institution_code"))
    col_coll   <- pick_col(c("collectionCode", "collection_code"))
    
    # Base extraction
    recordID <- if (!is.null(col_record)) as.character(df[[col_record]]) else NA_character_
    scientificName <- if (!is.null(col_sci)) as.character(df[[col_sci]]) else species_name
    eventDate <- if (!is.null(col_date)) as.character(df[[col_date]]) else NA_character_
    
    year <- if (!is.null(col_year)) {
      suppressWarnings(as.integer(df[[col_year]]))
    } else {
      suppressWarnings(as.integer(substr(as.character(eventDate), 1, 4)))
    }
    
    lat <- if (!is.null(col_lat)) suppressWarnings(as.numeric(df[[col_lat]])) else NA_real_
    lon <- if (!is.null(col_lon)) suppressWarnings(as.numeric(df[[col_lon]])) else NA_real_
    
    # Point fallback: parse "lat,lon" if decimal lat/lon were not provided or are missing.
    if (!is.null(col_point)) {
      x <- as.character(df[[col_point]])
      m <- stringr::str_match(x, "^\\s*([-+]?\\d+(?:\\.\\d+)?)\\s*,\\s*([-+]?\\d+(?:\\.\\d+)?)\\s*$")
      lat2 <- suppressWarnings(as.numeric(m[, 2]))
      lon2 <- suppressWarnings(as.numeric(m[, 3]))
      
      # If values look swapped (rare), swap them back
      swap <- !is.na(lat2) & !is.na(lon2) & (abs(lat2) > 90 & abs(lon2) <= 90)
      if (any(swap, na.rm = TRUE)) {
        tmp <- lat2[swap]
        lat2[swap] <- lon2[swap]
        lon2[swap] <- tmp
      }
      
      lat[is.na(lat)] <- lat2[is.na(lat)]
      lon[is.na(lon)] <- lon2[is.na(lon)]
    }
    
    license <- if (!is.null(col_lic)) as.character(df[[col_lic]]) else NA_character_
    cuim <- if (!is.null(col_cuim)) suppressWarnings(as.numeric(df[[col_cuim]])) else NA_real_
    cp <- if (!is.null(col_cp)) as.character(df[[col_cp]]) else NA_character_
    iv <- if (!is.null(col_iv)) as.character(df[[col_iv]]) else NA_character_
    idby <- if (!is.null(col_idby)) as.character(df[[col_idby]]) else NA_character_
    
    basis <- if (!is.null(col_basis)) as.character(df[[col_basis]]) else NA_character_
    rank <- if (!is.null(col_rank)) as.character(df[[col_rank]]) else NA_character_
    occst <- if (!is.null(col_occst)) as.character(df[[col_occst]]) else NA_character_
    dk <- if (!is.null(col_dk)) as.character(df[[col_dk]]) else NA_character_
    dn <- if (!is.null(col_dn)) as.character(df[[col_dn]]) else NA_character_
    pok <- if (!is.null(col_pok)) as.character(df[[col_pok]]) else NA_character_
    inst <- if (!is.null(col_inst)) as.character(df[[col_inst]]) else NA_character_
    coll <- if (!is.null(col_coll)) as.character(df[[col_coll]]) else NA_character_
    
    out <- data.frame(
      recordID = recordID,
      scientificName = scientificName,
      eventDate = eventDate,
      year = year,
      decimalLatitude = lat,
      decimalLongitude = lon,
      license = license,
      coordinateUncertaintyInMeters = cuim,
      coordinatePrecision = cp,
      identificationVerificationStatus = iv,
      identifiedBy = idby,
      basisOfRecord = basis,
      taxonRank = rank,
      occurrenceStatus = occst,
      datasetKey = dk,
      datasetName = dn,
      publishingOrgKey = pok,
      institutionCode = inst,
      collectionCode = coll,
      stringsAsFactors = FALSE
    )
    
    # Keep gridReference in the raw standardisation environment if present (useful for debugging),
    # but do not rely on it for coordinate conversion at this stage.
    if (!is.null(col_grid) && !"gridReference" %in% names(out)) {
      out$gridReference <- as.character(df[[col_grid]])
    }
    
    out
  }
  
  nbn_records_ws_total <- function(guid) {
    u <- paste0(
      "https://records-ws.nbnatlas.org/occurrences/search?",
      nbn_build_query(list(
        q = paste0("lsid:", guid),
        fq = '-occurrence_status:"absent"',
        pageSize = 0
      ))
    )
    raw <- tryCatch(jsonlite::fromJSON(u), error = function(e) e)
    if (inherits(raw, "error")) return(NA_integer_)
    suppressWarnings(as.integer(raw$totalRecords))
  }
  
  nbn_records_ws_download <- function(guid) {
    # Prefer the download endpoint for bulk export.
    # We set dwcHeaders=true for Darwin Core names and qa=none to avoid the extra assertions file.
    # Note: the NBN occurrence download web service is capped to 500,000 records per download.
    # For very common taxa, this may be a truncated export (flagged later via totalRecords).
    
    work_root <- file.path(get_checkpoint_root(repo_root), "nbn_work")
    dir.create(work_root, recursive = TRUE, showWarnings = FALSE)
    
    dl_base <- "https://records-ws.nbnatlas.org/occurrences/index/download"
    params <- list(
      q = paste0("lsid:", guid),
      fq = '-occurrence_status:"absent"',
      email = nbn_email,
      reasonTypeId = download_reason_id,
      fileType = "csv",
      dwcHeaders = "true",
      qa = "none"
    )
    
    dl_url <- paste0(dl_base, "?", nbn_build_query(params))
    
    ts <- format(Sys.time(), "%Y%m%d%H%M%S")
    zip_path <- file.path(work_root, paste0("nbn_download_", slug, "_", ts, ".zip"))
    unzip_dir <- file.path(work_root, paste0("nbn_download_", slug, "_", ts, "_unzipped"))
    
    # Longer downloads are common for large taxa; raise timeout locally around the download call.
    old_timeout <- getOption("timeout")
    options(timeout = max(as.integer(old_timeout), as.integer(nbn_download_timeout_s)))
    on.exit(options(timeout = old_timeout), add = TRUE)
    
    download_one <- function(dest) {
      # Prefer curl if available (more robust on Windows for large files).
      if (requireNamespace("curl", quietly = TRUE)) {
        res <- tryCatch({
          curl::curl_download(dl_url, destfile = dest, quiet = TRUE, mode = "wb")
          0L
        }, error = function(e) 1L)
        return(res)
      }
      
      rc <- tryCatch(
        utils::download.file(dl_url, destfile = dest, mode = "wb", quiet = TRUE, method = "libcurl"),
        warning = function(w) 1L,
        error = function(e) 1L
      )
      
      if (!identical(rc, 0L)) {
        rc2 <- tryCatch(
          utils::download.file(dl_url, destfile = dest, mode = "wb", quiet = TRUE),
          warning = function(w) 1L,
          error = function(e) 1L
        )
        rc <- rc2
      }
      
      rc
    }
    
    # Robust download handling:
    #   - handle timeouts / transient failures with one retry
    #   - treat non-zero return codes, missing files, or very small files as a failure
    rc <- download_one(zip_path)
    if (!identical(rc, 0L)) {
      unlink(zip_path, force = TRUE)
      Sys.sleep(2)
      rc <- download_one(zip_path)
    }
    
    if (!identical(rc, 0L) || !file.exists(zip_path) || is.na(file.info(zip_path)$size) || file.info(zip_path)$size < 200) {
      if (file.exists(zip_path)) unlink(zip_path, force = TRUE)
      stop("NBN records-ws download failed or returned an empty stub (rc=", rc, ").")
    }
    
    on.exit({
      unlink(zip_path, force = TRUE)
      unlink(unzip_dir, recursive = TRUE, force = TRUE)
    }, add = TRUE)
    
    if (nbn_is_zip_file(zip_path)) {
      dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
      utils::unzip(zip_path, exdir = unzip_dir)
      
      files <- list.files(unzip_dir, pattern = "\\.(csv|txt|tsv)$", full.names = TRUE,
                          ignore.case = TRUE, recursive = TRUE)
      
      f <- nbn_choose_occurrence_data_file(files)
      if (is.na(f) || !file.exists(f)) stop("NBN download unzip succeeded but could not find an occurrence CSV/TXT inside.")
      
      df <- nbn_read_download_file(f)
      return(nbn_standardise_ws_raw(df))
    }
    
    # Sometimes the endpoint returns CSV directly rather than a zip; handle that too.
    # However, the endpoint can also return an HTML error page; catch that early rather than parsing garbage.
    snip <- tryCatch(rawToChar(readBin(zip_path, "raw", n = 200)), error = function(e) "")
    if (nzchar(snip) && grepl("<html|service unavailable|request rejected|error", snip, ignore.case = TRUE)) {
      stop("NBN records-ws download returned a non-zip HTML response; treating as a failure.")
    }
    
    df <- nbn_read_download_file(zip_path)
    nbn_standardise_ws_raw(df)
  }
  
  nbn_records_ws_search_paged <- function(guid, page_size = 1000L) {
    # JSON paging fallback. Slower than downloads, and may not support deep paging for very large result sets.
    start <- 0L
    all <- list()
    
    repeat {
      u <- paste0(
        "https://records-ws.nbnatlas.org/occurrences/search?",
        nbn_build_query(list(
          q = paste0("lsid:", guid),
          fq = '-occurrence_status:"absent"',
          pageSize = page_size,
          startIndex = start
        ))
      )
      
      raw <- tryCatch(jsonlite::fromJSON(u), error = function(e) e)
      if (inherits(raw, "error")) stop("NBN records-ws search failed: ", conditionMessage(raw))
      
      occ <- raw$occurrences
      
      # If the service reports totalRecords but returns no occurrences before we reach it, treat this as an early stop.
      total <- suppressWarnings(as.integer(raw$totalRecords))
      if (is.null(occ) || length(occ) == 0) {
        if (!is.na(total) && isTRUE(total > start)) {
          stop("NBN records-ws paging ended early at startIndex=", start, " but totalRecords=", total)
        }
        break
      }
      
      df <- as.data.frame(occ, stringsAsFactors = FALSE)
      all[[length(all) + 1]] <- df
      
      got <- start + nrow(df)
      message("[NBN] records-ws paging: ", got, " / ", raw$totalRecords)
      
      if (!is.numeric(raw$totalRecords) || got >= raw$totalRecords) break
      
      start <- got
      Sys.sleep(pause_s)
    }
    
    nbn_standardise_ws_raw(dplyr::bind_rows(all))
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
  #   - For non-empty cached files: trust the cache.
  #   - For empty cached files: trust the cache only if nbn_state says complete=TRUE.
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
    } else if (nrow(nbn_cached) > 0) {
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
      
      attr(nbn_clean, "nbn_status") <- list(
        state = "complete",
        note = if (!is.null(nbn_state$note)) nbn_state$note else NA_character_,
        totalRecords = if (!is.null(nbn_state$totalRecords)) nbn_state$totalRecords else NA_integer_
      )
      
      return(nbn_clean)
      
    } else {
      # Empty cached file: only treat as final if we previously recorded this as a complete pull.
      if (isTRUE(nbn_state$complete)) {
        message("Found existing NBN clean file (EMPTY) and NBN state is complete; keeping: ", nbn_outfile)
        nbn_clean <- empty_nbn_clean()
        attr(nbn_clean, "nbn_status") <- list(state = "complete", note = nbn_state$note, totalRecords = nbn_state$totalRecords)
        return(nbn_clean)
      } else {
        message("Found existing NBN clean file (EMPTY) but NBN state is not complete; retrying NBN pull: ", nbn_outfile)
        use_cached <- FALSE
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # NBN taxon guard
  #   We proceed only if NBN taxonomy contains an exact, species-rank match.
  # ---------------------------------------------------------------------------
  nbn_use_ws <- FALSE
  nbn_guid <- NA_character_
  
  nbn_taxa <- tryCatch(search_taxa(species_name), error = function(e) e)
  
  if (inherits(nbn_taxa, "error")) {
    msg <- conditionMessage(nbn_taxa)
    message(
      "[NBN] galah taxon lookup failed (will use NBN web services directly).\n",
      "      Error: ", msg
    )
    nbn_use_ws <- TRUE
  } else {
    message("NBN taxon search (top hit):")
    
    tryCatch(
      {
        ensure_safe_na_print()
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
      },
      error = function(e) {
        ensure_safe_na_print()
        message("[NBN] NOTE: could not print taxon search preview (", conditionMessage(e), "). Continuing.")
      }
    )
    
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
        
        nbn_state$complete <- TRUE
        nbn_state$last_updated <- as.character(Sys.time())
        nbn_state$note <- "no_exact_species_match"
        nbn_state$last_error <- NA_character_
        nbn_state$totalRecords <- 0L
        safe_saveRDS(nbn_state, nbn_state_file)
        
        attr(nbn_clean, "nbn_status") <- list(state = "complete", note = nbn_state$note, totalRecords = nbn_state$totalRecords)
        return(nbn_clean)
      }
      
      # If the taxonomy table includes a concept ID/guid, keep it for the records-ws fallback.
      id_col <- intersect(c("taxon_concept_id", "taxonConceptId", "guid"), names(nbn_exact))
      if (length(id_col) > 0) {
        nbn_guid <- as.character(nbn_exact[[id_col[1]]][1])
      }
      
    } else {
      message("[NBN] Taxon table format unexpected; proceeding to attempt pull.")
    }
  }
  
  max_retries <- 5
  retry_base_wait_s <- 10
  
  nbn_raw <- NULL
  last_err <- NA_character_
  
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
  
  if (!isTRUE(nbn_use_ws)) {
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
      
      last_err <- conditionMessage(nbn_raw_try)
      
      # Deterministic taxonomy parsing bug (non-rectangular taxa response); switch to web services.
      if (grepl("Can't recycle `id`|vctrs::data_frame\\(|synonymComplete", last_err, ignore.case = TRUE)) {
        message(
          "[NBN] galah occurrence pull hit a deterministic taxonomy parsing error.\n",
          "      Switching to NBN web services for this species.\n",
          "      Error: ", last_err
        )
        nbn_use_ws <- TRUE
        break
      }
      
      message(
        "NBN combined (core+QA) pull failed (attempt ", attempt, "/", max_retries, "): ",
        last_err,
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
        last_err <- conditionMessage(core_try)
        
        if (grepl("Can't recycle `id`|vctrs::data_frame\\(|synonymComplete", last_err, ignore.case = TRUE)) {
          message(
            "[NBN] galah occurrence pull hit a deterministic taxonomy parsing error.\n",
            "      Switching to NBN web services for this species.\n",
            "      Error: ", last_err
          )
          nbn_use_ws <- TRUE
          break
        }
        
        wait_s <- retry_base_wait_s * attempt
        message(
          "NBN core pull also failed (attempt ", attempt, "/", max_retries, "): ",
          last_err,
          " | waiting ", wait_s, "s then retrying..."
        )
        Sys.sleep(wait_s)
      }
    }
  }
  
  if (isTRUE(nbn_use_ws) && is.null(nbn_raw)) {
    
    message("[NBN] Using NBN web services fallback for: ", species_name)
    
    if (is.na(nbn_guid) || !nzchar(nbn_guid)) {
      res <- tryCatch(nbn_species_ws_search(species_name), error = function(e) e)
      
      if (inherits(res, "error")) {
        nbn_state$complete <- FALSE
        nbn_state$last_updated <- as.character(Sys.time())
        nbn_state$note <- "taxon_lookup_failed_species_ws"
        nbn_state$last_error <- conditionMessage(res)
        safe_saveRDS(nbn_state, nbn_state_file)
        stop("NBN species-ws lookup failed for: ", species_name, " | ", conditionMessage(res))
      }
      
      nbn_guid <- nbn_pick_guid(res, species_name)
      
      if (is.na(nbn_guid) || !nzchar(nbn_guid)) {
        message(
          "[NBN] No exact species-rank match for '", species_name, "' in species-ws.\n",
          "      Skipping NBN pull and writing an empty output so the pipeline can continue."
        )
        
        nbn_clean <- empty_nbn_clean()
        readr::write_csv(nbn_clean, nbn_outfile)
        message("Saved NBN clean file (EMPTY): ", nbn_outfile)
        
        nbn_state$complete <- TRUE
        nbn_state$last_updated <- as.character(Sys.time())
        nbn_state$note <- "no_exact_species_match_species_ws"
        nbn_state$last_error <- NA_character_
        nbn_state$totalRecords <- 0L
        safe_saveRDS(nbn_state, nbn_state_file)
        
        attr(nbn_clean, "nbn_status") <- list(state = "complete", note = nbn_state$note, totalRecords = nbn_state$totalRecords)
        return(nbn_clean)
      }
      
      message("[NBN] species-ws match GUID: ", nbn_guid)
    }
    
    # Record totalRecords for QA (also used to flag likely-truncated downloads)
    nbn_state$guid <- nbn_guid
    nbn_state$totalRecords <- nbn_records_ws_total(nbn_guid)
    if (!is.na(nbn_state$totalRecords) && nbn_state$totalRecords > 500000L) {
      message(
        "\n[NBN] NOTE: totalRecords=", nbn_state$totalRecords,
        " for ", species_name, ". The records-ws download endpoint is capped to 500,000 rows per download.\n",
        "If you need full coverage for this taxon, it must be retrieved in multiple filtered downloads (e.g. by year ranges).\n"
      )
    }
    
    # Try bulk download first; if it errors, fall back to paged JSON search.
    nbn_raw <- tryCatch(nbn_records_ws_download(nbn_guid), error = function(e) e)
    
    if (inherits(nbn_raw, "error")) {
      message("[NBN] records-ws download failed; trying JSON paging fallback. Error: ", conditionMessage(nbn_raw))
      nbn_raw <- tryCatch(nbn_records_ws_search_paged(nbn_guid), error = function(e) e)
    }
    
    if (inherits(nbn_raw, "error")) {
      nbn_state$complete <- FALSE
      nbn_state$last_updated <- as.character(Sys.time())
      nbn_state$note <- "records_ws_failed"
      nbn_state$last_error <- conditionMessage(nbn_raw)
      safe_saveRDS(nbn_state, nbn_state_file)
      stop("NBN records-ws fallback failed for: ", species_name, " | ", conditionMessage(nbn_raw))
    }
  }
  
  if (is.null(nbn_raw)) {
    nbn_state$complete <- FALSE
    nbn_state$last_updated <- as.character(Sys.time())
    nbn_state$note <- "download_failed_after_retries"
    nbn_state$last_error <- if (!is.na(last_err)) last_err else "NBN pull failed after retries"
    safe_saveRDS(nbn_state, nbn_state_file)
    stop("NBN pull failed after retries for: ", species_name)
  }
  
  message("NBN raw rows: ", nrow(nbn_raw))
  
  lic_col <- if ("dcterms:license" %in% names(nbn_raw)) "dcterms:license" else if ("license" %in% names(nbn_raw)) "license" else "license"
  
  if (!"coordinateUncertaintyInMeters" %in% names(nbn_raw)) nbn_raw$coordinateUncertaintyInMeters <- NA_real_
  if (!"coordinatePrecision" %in% names(nbn_raw))          nbn_raw$coordinatePrecision <- NA_character_
  if (!"identificationVerificationStatus" %in% names(nbn_raw)) nbn_raw$identificationVerificationStatus <- NA_character_
  if (!"identifiedBy" %in% names(nbn_raw))                 nbn_raw$identifiedBy <- NA_character_
  
  # Schema alignment / provenance fields (often NA for NBN)
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
    
    tryCatch(
      {
        ensure_safe_na_print()
        nbn_raw %>%
          count(.data[[lic_col]], sort = TRUE) %>%
          mutate(prop = n / sum(n)) %>%
          print(n = 10)
      },
      error = function(e) {
        ensure_safe_na_print()
        message("[NBN] NOTE: could not print licence breakdown (", conditionMessage(e), "). Continuing.")
      }
    )
  } else {
    message("[NBN] No rows returned (0 UK records is plausible for non-native taxa).")
  }
  
  # Harmonise name variants that show up between galah and records-ws exports
  if (!"decimalLongitude" %in% names(nbn_raw) && "decimal_longitude" %in% names(nbn_raw)) nbn_raw$decimalLongitude <- nbn_raw$decimal_longitude
  if (!"decimalLatitude"  %in% names(nbn_raw) && "decimal_latitude"  %in% names(nbn_raw)) nbn_raw$decimalLatitude  <- nbn_raw$decimal_latitude
  
  nbn_clean <- nbn_raw %>%
    transmute(
      source = "NBN",
      species = species_name,
      recordID = as.character(recordID),
      lon = suppressWarnings(as.numeric(decimalLongitude)),
      lat = suppressWarnings(as.numeric(decimalLatitude)),
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
  
  nbn_state$complete <- TRUE
  nbn_state$last_updated <- as.character(Sys.time())
  
  # Flag likely-truncated downloads for very common taxa (the endpoint is capped to 500k).
  # We keep the file (it is still useful), but record the condition for a later QA step.
  if (!is.na(nbn_state$totalRecords) && nbn_state$totalRecords > 500000L) {
    nbn_state$note <- "complete_totalRecords_gt_500k_possible_truncation"
  } else if (nrow(nbn_clean) == 0) {
    nbn_state$note <- "complete_zero_records"
  } else {
    nbn_state$note <- "complete"
  }
  
  if (isTRUE(nbn_use_ws)) nbn_state$note <- paste0(nbn_state$note, "_via_records_ws")
  
  nbn_state$last_error <- NA_character_
  safe_saveRDS(nbn_state, nbn_state_file)
  
  attr(nbn_clean, "nbn_status") <- list(
    state = "complete",
    note = nbn_state$note,
    totalRecords = nbn_state$totalRecords,
    guid = nbn_state$guid
  )
  
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
    
    tryCatch(
      {
        ensure_safe_na_print()
        print(dplyr::slice_head(unexpected_tbl, n = 10), n = 10)
      },
      error = function(e) {
        ensure_safe_na_print()
        message("[NBN] NOTE: could not print unexpected licence table (", conditionMessage(e), "). Continuing.")
      }
    )
    
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
  
  tryCatch(
    {
      ensure_safe_na_print()
      nbn_clean %>% count(licence, sort = TRUE) %>% print(n = 10)
    },
    error = function(e) {
      ensure_safe_na_print()
      message("[NBN] NOTE: could not print final licence table (", conditionMessage(e), "). Continuing.")
    }
  )
  
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
                                 gbif_download_on_search_error = TRUE,
                                 skip_species_if_complete = TRUE,
                                 cleanup_gbif_work_files = TRUE,
                                 nbn_download_timeout_s = 3600L) {
  
  if (missing(nbn_email) || is.null(nbn_email) || !nzchar(nbn_email)) {
    stop("Please provide nbn_email (the email associated with your NBN Atlas account).")
  }
  
  repo_root <- get_repo_root()
  group_dir2 <- normalise_group_dir(group_dir)
  
  species_complete <- function(sp) {
    slug <- slugify_species(sp)
    
    gbif_out_root <- if (nzchar(group_dir2)) file.path(repo_root, "data", "raw", "gbif", group_dir2) else file.path(repo_root, "data", "raw", "gbif")
    nbn_out_root  <- if (nzchar(group_dir2)) file.path(repo_root, "data", "raw", "nbn",  group_dir2) else file.path(repo_root, "data", "raw", "nbn")
    
    gbif_out_dir <- if (isTRUE(species_subdir)) file.path(gbif_out_root, slug) else gbif_out_root
    nbn_out_dir  <- if (isTRUE(species_subdir)) file.path(nbn_out_root,  slug) else nbn_out_root
    
    gbif_csv <- file.path(gbif_out_dir, paste0("gbif_", slug, "_clean.csv"))
    nbn_csv  <- file.path(nbn_out_dir,  paste0("nbn_",  slug, "_clean.csv"))
    
    ckpt_root <- get_checkpoint_root(repo_root)
    gbif_ckpt <- file.path(ckpt_root, "gbif", paste0("gbif_pull_checkpoint_", slug, ".rds"))
    nbn_state <- file.path(ckpt_root, "nbn",  paste0("nbn_state_", slug, ".rds"))
    
    gbif_ok <- FALSE
    if (file.exists(gbif_csv) && file.exists(gbif_ckpt)) {
      x <- tryCatch(readRDS(gbif_ckpt), error = function(e) NULL)
      gbif_ok <- is.list(x) && isTRUE(x$complete)
    }
    
    nbn_ok <- FALSE
    if (file.exists(nbn_csv) && file.exists(nbn_state)) {
      y <- tryCatch(readRDS(nbn_state), error = function(e) NULL)
      nbn_ok <- is.list(y) && isTRUE(y$complete)
    }
    
    isTRUE(gbif_ok && nbn_ok)
  }
  
  out <- vector("list", length(species_names))
  names(out) <- species_names
  
  gbif_incomplete <- list()
  nbn_incomplete  <- list()
  
  for (i in seq_along(species_names)) {
    sp <- species_names[i]
    
    if (isTRUE(skip_species_if_complete) && isTRUE(species_complete(sp))) {
      message("\n============================================================")
      message("Pulling raw occurrence outputs for: ", sp)
      message("============================================================\n")
      message("[SKIP] Outputs already complete for GBIF + NBN; moving on.")
      next
    }
    
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
      gbif_download_on_search_error = gbif_download_on_search_error,
      cleanup_gbif_work_files = cleanup_gbif_work_files
    )
    
    st <- attr(gbif_clean, "gbif_status")
    if (!is.null(st) && !identical(st$state, "complete")) {
      gbif_incomplete[[sp]] <- st
    }
    
    nbn_clean <- tryCatch(
      pull_nbn_clean(
        species_name = sp,
        group_dir = group_dir,
        nbn_email = nbn_email,
        download_reason_id = download_reason_id,
        expected_licences_nbn = expected_licences_nbn,
        use_cache = use_cache,
        pause_s = pause_s,
        species_subdir = species_subdir,
        nbn_download_timeout_s = nbn_download_timeout_s
      ),
      error = function(e) {
        msg <- conditionMessage(e)
        
        reason <- if (grepl("OAuth error|authentication required|HTTP 403", msg, ignore.case = TRUE)) {
          "authentication_required"
        } else {
          "error"
        }
        
        message(
          "\n[NBN][INCOMPLETE] Failed to pull NBN records for ", sp, ".\n",
          "  Error: ", msg, "\n",
          "Continuing to the next species.\n"
        )
        
        if (identical(reason, "authentication_required")) {
          message(
            "[NBN] This looks like an authentication problem. In an interactive session, run:\n",
            "  library(galah)\n",
            "  galah_config(atlas = \"United Kingdom\", email = \"", nbn_email, "\", verbose = FALSE)\n",
            "  galah_login()\n"
          )
        }
        
        # Write an empty output so downstream steps don't break.
        # However, if a non-empty clean file already exists (e.g. the pull succeeded and only a final print failed),
        # do not overwrite it with an empty file.
        group_dir3 <- normalise_group_dir(group_dir)
        slug <- slugify_species(sp)
        nbn_out_root <- if (nzchar(group_dir3)) file.path(repo_root, "data", "raw", "nbn", group_dir3) else file.path(repo_root, "data", "raw", "nbn")
        nbn_out_dir <- if (isTRUE(species_subdir)) file.path(nbn_out_root, slug) else nbn_out_root
        dir.create(nbn_out_dir, recursive = TRUE, showWarnings = FALSE)
        nbn_outfile <- file.path(nbn_out_dir, paste0("nbn_", slug, "_clean.csv"))
        
        existing_nonempty <- FALSE
        existing_tbl <- NULL
        
        if (file.exists(nbn_outfile)) {
          existing_tbl <- tryCatch(readr::read_csv(nbn_outfile, show_col_types = FALSE), error = function(e2) NULL)
          if (!is.null(existing_tbl) && is.data.frame(existing_tbl) && nrow(existing_tbl) > 0) {
            existing_nonempty <- TRUE
          }
        }
        
        if (isTRUE(existing_nonempty)) {
          message("Found existing NBN clean file written before the error; keeping: ", nbn_outfile)
          
          # Mark NBN state as complete so the existing file is treated as final on the next run.
          ckpt_root <- get_checkpoint_root(repo_root)
          nbn_state_file <- file.path(ckpt_root, "nbn", paste0("nbn_state_", slug, ".rds"))
          dir.create(dirname(nbn_state_file), recursive = TRUE, showWarnings = FALSE)
          nbn_state <- list(
            complete = TRUE,
            last_updated = as.character(Sys.time()),
            note = "complete_existing_csv_kept_after_error",
            last_error = msg,
            totalRecords = NA_integer_,
            guid = NA_character_
          )
          safe_saveRDS(nbn_state, nbn_state_file)
          
          # Return the existing table, with an informative status attribute for the wrapper summaries.
          attr(existing_tbl, "nbn_status") <- list(
            state = "complete",
            reason = "kept_existing_nonempty_csv_after_error",
            error = msg
          )
          
          return(existing_tbl)
        }
        
        nbn_clean_fallback <- tibble::tibble(
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
          basisOfRecord = character(),
          taxonRank = character(),
          occurrenceStatus = character(),
          datasetKey = character(),
          datasetName = character(),
          publishingOrgKey = character(),
          institutionCode = character(),
          collectionCode = character()
        )
        
        readr::write_csv(nbn_clean_fallback, nbn_outfile)
        message("Saved NBN clean file (EMPTY): ", nbn_outfile)
        
        # Mark NBN state as incomplete so an empty file does not look "done" on the next run
        ckpt_root <- get_checkpoint_root(repo_root)
        nbn_state_file <- file.path(ckpt_root, "nbn", paste0("nbn_state_", slug, ".rds"))
        dir.create(dirname(nbn_state_file), recursive = TRUE, showWarnings = FALSE)
        nbn_state <- list(
          complete = FALSE,
          last_updated = as.character(Sys.time()),
          note = "incomplete",
          last_error = msg,
          totalRecords = NA_integer_,
          guid = NA_character_
        )
        safe_saveRDS(nbn_state, nbn_state_file)
        
        attr(nbn_clean_fallback, "nbn_status") <- list(
          state = "incomplete",
          reason = reason,
          error = msg
        )
        
        nbn_clean_fallback
      }
    )
    
    nst <- attr(nbn_clean, "nbn_status")
    if (!is.null(nst) && !identical(nst$state, "complete")) {
      nbn_incomplete[[sp]] <- nst
    }
    
    out[[i]] <- list(gbif_clean = gbif_clean, nbn_clean = nbn_clean)
  }
  
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
  
  if (length(nbn_incomplete) > 0) {
    message("\n==================== NBN WARNING ====================")
    message("Some NBN pulls did not complete (often authentication or temporary service issues).")
    for (nm in names(nbn_incomplete)) {
      s <- nbn_incomplete[[nm]]
      line <- paste0(" - ", nm, ": ", s$reason)
      if (!is.null(s$error) && !is.na(s$error) && nzchar(s$error)) line <- paste0(line, " | ", s$error)
      message(line)
    }
    message("====================================================\n")
  }
  
  invisible(out)
}
