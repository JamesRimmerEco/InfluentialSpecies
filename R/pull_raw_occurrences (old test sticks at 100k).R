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
#   Save outputs for each species to species subfolders in raw:
#     data/raw/gbif/<group_dir>/<slug>/gbif_<slug>_clean.csv
#     data/raw/nbn/<group_dir>/<slug>/nbn_<slug>_clean.csv
#
#   Save checkpoints to:
#     data/_checkpoints/gbif/gbif_pull_checkpoint_<slug>.rds
#     data/_checkpoints/nbn/nbn_pull_checkpoint_<slug>.rds
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

# ---- Helper: slugify a species name ------------------------------------------
slugify_species <- function(species_name) {
  slug <- str_replace_all(tolower(species_name), "[^a-z0-9]+", "_")
  slug <- str_replace_all(slug, "^_+|_+$", "")
  slug
}

# "Slugify" - just means creating a user-friendly descriptive domain, easier to read
# than codes.

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
                            group_dir = "wasps",
                            pause_s = 0.25,
                            page_size = 1000,
                            max_records = Inf,
                            expected_licences_gbif = c("CC0_1_0", "CC_BY_4_0", "CC_BY_NC_4_0"),
                            use_cache = TRUE,
                            species_subdir = TRUE) {
  
  repo_root <- get_repo_root()
  slug <- slugify_species(species_name)
  
  # Output dirs
  gbif_out_root <- file.path(repo_root, "data", "raw", "gbif", group_dir)
  ckpt_root     <- file.path(repo_root, "data", "_checkpoints")
  gbif_ckpt_dir <- file.path(ckpt_root, "gbif")
  
  dir.create(gbif_out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(gbif_ckpt_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Species subfolder (NEW)
  gbif_out_dir <- if (isTRUE(species_subdir)) file.path(gbif_out_root, slug) else gbif_out_root
  dir.create(gbif_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Clean output file (used as a cache when re-running)
  gbif_outfile <- file.path(gbif_out_dir, paste0("gbif_", slug, "_clean.csv"))
  
  # Licence normaliser (GBIF uses URLs with 'licenses' in the path)
  lic_normalise_gbif <- function(x) {
    dplyr::case_when(
      stringr::str_detect(x, "publicdomain/zero/1.0") ~ "CC0_1_0",
      stringr::str_detect(x, "licenses/by/4.0") ~ "CC_BY_4_0",
      stringr::str_detect(x, "licenses/by-nc/4.0") ~ "CC_BY_NC_4_0",
      TRUE ~ NA_character_
    )
  }
  
  # Columns that MUST exist in cached outputs to be considered "current schema"
  required_cache_cols <- c(
    "coordinateUncertaintyInMeters",
    "identificationVerificationStatus",
    "issues"
  )
  
  # ---------------------------------------------------------------------------
  # GBIF taxon resolution
  # ---------------------------------------------------------------------------
  bb <- name_backbone(name = species_name)
  taxon_key <- bb$usageKey
  message("GBIF match: ", bb$scientificName, " (usageKey=", taxon_key, ", matchType=", bb$matchType, ")")
  
  # ---------------------------------------------------------------------------
  # GBIF occurrence pull (cached if available + schema matches)
  # ---------------------------------------------------------------------------
  use_cached <- isTRUE(use_cache) && file.exists(gbif_outfile)
  
  if (use_cached) {
    gbif_cached <- readr::read_csv(gbif_outfile, show_col_types = FALSE)
    missing_cols <- setdiff(required_cache_cols, names(gbif_cached))
    
    if (length(missing_cols) > 0) {
      message(
        "Found existing GBIF clean file but it is missing new QA columns: ",
        paste(missing_cols, collapse = ", "),
        "\nRe-pulling from GBIF: ", gbif_outfile
      )
      use_cached <- FALSE
    } else {
      message("Found existing GBIF clean file; skipping API pull: ", gbif_outfile)
      
      gbif_clean <- gbif_cached %>%
        mutate(
          licence = as.character(licence),
          licence_raw = as.character(licence_raw),
          licence_expected = !is.na(licence) & licence %in% expected_licences_gbif,
          identificationVerificationStatus = as.character(identificationVerificationStatus),
          issues = as.character(issues)
        )
    }
  }
  
  if (!use_cached) {
    
    ckpt_file <- file.path(gbif_ckpt_dir, paste0("gbif_pull_checkpoint_", slug, ".rds"))
    
    max_retries <- 5
    retry_base_wait_s <- 10
    
    start <- 0
    all_pages <- list()
    total_expected <- NA_integer_
    
    # If a checkpoint exists, resume instead of starting from scratch
    if (file.exists(ckpt_file)) {
      ckpt <- readRDS(ckpt_file)
      start <- ckpt$start
      all_pages <- ckpt$all_pages
      total_expected <- ckpt$total_expected
      
      message("Resuming GBIF pull from start = ", start,
              " (pages already stored: ", length(all_pages), ").")
    }
    
    repeat {
      Sys.sleep(pause_s)
      
      res <- NULL
      for (attempt in seq_len(max_retries)) {
        res <- tryCatch(
          occ_search(
            taxonKey = taxon_key,
            continent = region_scope,
            hasCoordinate = TRUE,
            limit = page_size,
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
        saveRDS(list(start = start, all_pages = all_pages, total_expected = total_expected), ckpt_file)
        stop(res)
      }
      
      if (is.na(total_expected)) total_expected <- res$meta$count
      if (length(res$data) == 0) break
      
      all_pages[[length(all_pages) + 1]] <- res$data
      
      saveRDS(list(start = start, all_pages = all_pages, total_expected = total_expected), ckpt_file)
      
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
    # Basic screening + essential fields (UPDATED: include QA/certainty fields)
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
  }
  
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
  
  return(gbif_clean)
}

# ==============================================================================
# NBN pull (UK Atlas) ----------------------------------------------------------
# ==============================================================================

pull_nbn_clean <- function(species_name,
                           group_dir = "wasps",
                           nbn_email,
                           download_reason_id = 17,
                           expected_licences_nbn = c("OGL", "CC0", "CC-BY", "CC-BY-NC"),
                           use_cache = TRUE,
                           pause_s = 0.25,
                           species_subdir = TRUE) {
  
  repo_root <- get_repo_root()
  slug <- slugify_species(species_name)
  
  # Output dirs
  nbn_out_root <- file.path(repo_root, "data", "raw", "nbn", group_dir)
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
    "identifiedBy"
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
        "Found existing NBN clean file but it is missing new QA columns: ",
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
          coordinatePrecision = as.character(coordinatePrecision)
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
        identifiedBy = as.character(identifiedBy)
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
                                 group_dir = "wasps",
                                 region_scope = "EUROPE",
                                 pause_s = 0.25,
                                 page_size = 1000,
                                 max_records = Inf,
                                 nbn_email,
                                 download_reason_id = 17,
                                 expected_licences_gbif = c("CC0_1_0", "CC_BY_4_0", "CC_BY_NC_4_0"),
                                 expected_licences_nbn  = c("OGL", "CC0", "CC-BY", "CC-BY-NC"),
                                 use_cache = TRUE,
                                 species_subdir = TRUE) {
  
  if (missing(nbn_email) || is.null(nbn_email) || !nzchar(nbn_email)) {
    stop("Please provide nbn_email (the email associated with your NBN Atlas account).")
  }
  
  out <- vector("list", length(species_names))
  names(out) <- species_names
  
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
      species_subdir = species_subdir
    )
    
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
  
  invisible(out)
}
