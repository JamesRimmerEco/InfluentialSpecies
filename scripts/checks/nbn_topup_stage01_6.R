# scripts/checks/nbn_topup_stage01_6.R
#
# InfluentialSpecies — Stage 1.6 (NBN top-up / repair)
#
# Purpose
#   Stage 1.5 tells you which species are missing / risky / cap-risk after a Stage 1 raw pull.
#   Stage 1.6 takes the Stage 1.5 manifest and *repairs NBN only* by re-downloading NBN
#   occurrences in year-range chunks to avoid the ~500k cap, then appends + de-duplicates
#   into the existing per-species NBN clean CSV (same schema as Stage 1).
#
# Key design choices (to keep outputs from “exploding”)
#   - Updates only: data/raw/nbn/<group_dir>/nbn_<slug>_clean.csv
#   - Writes one run summary CSV: data/_meta/audits/nbn_topup_summary_<group_dir>_<timestamp>.csv
#   - Keeps resumable checkpoints under INFLUENTIAL_CHECKPOINT_ROOT (or a sensible default)
#   - No per-chunk CSV outputs unless you explicitly add them later
#
# What “auto-work” means here
#   For each species flagged in the Stage 1.5 manifest:
#     1) Resolve GUID via NBN species-ws
#     2) Estimate totalRecords via records-ws for a year range
#     3) If near cap, recursively split year range until each chunk is “safe” (or single-year)
#     4) Download each chunk via records-ws download endpoint (zip -> csv)
#     5) Standardise columns to the Stage 1 NBN schema
#     6) Append + dedupe on recordID
#     7) Write updated per-species clean CSV and checkpoint progress
#
# Notes / assumptions
#   - Dedupe key for NBN is recordID (assumed present; script fails fast if missing after standardisation).
#   - Chunking is done by numeric year ranges (Solr fq syntax: year:[YYYY TO YYYY]).
#     This avoids needing month/day filters.
#   - If a single year is still near-cap, the script will still download (likely truncated)
#     and will flag that chunk as “too_big_single_year”. Handling that would require a deeper split
#     by date or another facet (future enhancement).
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(jsonlite)
})

# ---- User options -------------------------------------------------------------

# The run folder you want to top up (must match Stage 1 outputs).
group_dir <- "home_run_true_list"

# Start year for “all years” top-up.
start_year <- 1900L
end_year   <- as.integer(format(Sys.Date(), "%Y"))

# NBN “cap” heuristics (used for splitting decisions)
nbn_cap_n          <- 500000L
nbn_near_cap_prop  <- 0.98  # split if totalRecords >= cap*0.98
max_split_depth    <- 20L   # guardrail

# Which species to process from the manifest:
# - topup_recommended == TRUE OR nbn_status indicates a hard problem
process_if <- function(manifest_row) {
  isTRUE(manifest_row$topup_recommended) ||
    (manifest_row$nbn_status %in% c("missing_file", "unreadable", "missing_required_columns"))
}

# NBN download identity (required by records-ws download endpoint)
nbn_email <- "jamesrimmer92@mail.com"

# NBN download reasonTypeId:
# - Your Stage 1 engine likely sets this; here we keep it explicit and easy to edit.
# - If you already know the correct ID for your account/org, set it here.
# - If you don’t, leave as-is and the service may still accept the request depending on config.
download_reason_id <- 10

# Optional: only do a dry run (plan chunks + print, no downloads)
dry_run <- FALSE

# Download robustness
dl_max_tries         <- 4L     # total attempts per chunk
dl_backoff_base_sec  <- 4L     # base sleep; grows with attempt
dl_min_zip_bytes     <- 1000L  # treat zips smaller than this as suspicious
dl_debug_keep_failed <- TRUE   # keep failed payloads (html/txt) for inspection

# ---- Repo paths ---------------------------------------------------------------

repo_root <- getwd()

meta_dir   <- file.path(repo_root, "data", "_meta")
audits_dir <- file.path(meta_dir, "audits")

nbn_run_dir <- file.path(repo_root, "data", "raw", "nbn", group_dir)
dir.create(nbn_run_dir, recursive = TRUE, showWarnings = FALSE)

timestamp_tag <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
out_summary_csv <- file.path(audits_dir, paste0("nbn_topup_summary_", group_dir, "_", timestamp_tag, ".csv"))

# ---- Checkpoint root ----------------------------------------------------------

get_checkpoint_root <- function() {
  x <- Sys.getenv("INFLUENTIAL_CHECKPOINT_ROOT")
  if (nzchar(x)) return(x)
  # sensible fallback
  file.path(repo_root, "data", "_checkpoints")
}

ckpt_root <- get_checkpoint_root()
dir.create(ckpt_root, recursive = TRUE, showWarnings = FALSE)

nbn_topup_ckpt_dir <- file.path(ckpt_root, "nbn_topup")
dir.create(nbn_topup_ckpt_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers -----------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

slugify_species <- function(species_name) {
  slug <- str_replace_all(tolower(species_name), "[^a-z0-9]+", "_")
  slug <- str_replace_all(slug, "^_+|_+$", "")
  slug
}

nbn_build_query <- function(params) {
  # Build a query string from a list, allowing multiple fq values.
  # params$fq can be a single string or a character vector.
  parts <- character(0)
  for (nm in names(params)) {
    v <- params[[nm]]
    if (is.null(v)) next
    if (length(v) == 0) next
    
    if (nm == "fq" && length(v) > 1) {
      for (vv in v) {
        parts <- c(parts, paste0("fq=", utils::URLencode(vv, reserved = TRUE)))
      }
    } else {
      parts <- c(parts, paste0(nm, "=", utils::URLencode(as.character(v[1]), reserved = TRUE)))
    }
  }
  paste(parts, collapse = "&")
}

http_head_status <- function(url) {
  if (!requireNamespace("curl", quietly = TRUE)) return(NA_integer_)
  h <- curl::new_handle(nobody = TRUE, followlocation = TRUE)
  curl::handle_setopt(h, connecttimeout = 30, timeout = 60)
  res <- tryCatch(curl::curl_fetch_memory(url, handle = h), error = function(e) NULL)
  if (is.null(res)) return(NA_integer_)
  suppressWarnings(as.integer(res$status_code))
}

looks_like_html <- function(path, n = 6L) {
  if (!file.exists(path)) return(FALSE)
  x <- tryCatch(readLines(path, n = n, warn = FALSE), error = function(e) character(0))
  if (length(x) == 0) return(FALSE)
  any(grepl("<!DOCTYPE html|<html|<head|<body", x, ignore.case = TRUE))
}

save_failed_payload_preview <- function(path, out_path) {
  if (!file.exists(path)) return(FALSE)
  ok <- tryCatch({ file.copy(path, out_path, overwrite = TRUE); TRUE }, error = function(e) FALSE)
  if (!isTRUE(ok)) return(FALSE)
  prev <- tryCatch(readLines(out_path, n = 80, warn = FALSE), error = function(e) character(0))
  if (length(prev)) {
    tryCatch(writeLines(prev, paste0(out_path, ".preview.txt")), error = function(e) NULL)
  }
  TRUE
}

# ---- Stage 1.5 manifest discovery --------------------------------------------

find_latest_manifest <- function(group_dir) {
  if (!dir.exists(audits_dir)) return(NA_character_)
  pat <- paste0("^audit_raw_pull_manifest_", stringr::fixed(group_dir), "_.+\\.csv$")
  hits <- list.files(audits_dir, pattern = pat, full.names = TRUE)
  if (length(hits) == 0) return(NA_character_)
  hits[which.max(file.info(hits)$mtime)]
}

manifest_path <- find_latest_manifest(group_dir)
if (is.na(manifest_path)) {
  stop(
    "Could not find a Stage 1.5 manifest for group_dir='", group_dir, "' in ",
    normalizePath(audits_dir, winslash = "/", mustWork = FALSE), "\n",
    "Run scripts/checks/audit_raw_pull_v1.R first (with do_nbn_online_probe=TRUE)."
  )
}

message("[1.6] Using manifest: ", normalizePath(manifest_path, winslash = "/", mustWork = FALSE))
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)

need_cols <- c("species", "slug", "nbn_status", "topup_recommended")
missing_cols <- setdiff(need_cols, names(manifest))
if (length(missing_cols) > 0) {
  stop("Manifest missing required columns: ", paste(missing_cols, collapse = ", "))
}

worklist <- manifest %>%
  mutate(
    .do = (topup_recommended %in% TRUE) |
      (nbn_status %in% c("missing_file", "unreadable", "missing_required_columns"))
  ) %>%
  filter(.do) %>%
  select(species, slug, nbn_status, topup_recommended)

message("[1.6] Worklist size: ", nrow(worklist), " species")

if (nrow(worklist) == 0) {
  message("[1.6] Nothing to do; exiting.")
  quit(save = "no")
}

# ---- NBN species-ws: resolve GUID (LSID) -------------------------------------

nbn_species_ws_search <- function(sp) {
  u <- paste0(
    "https://species-ws.nbnatlas.org/search?q=",
    utils::URLencode(sp, reserved = TRUE),
    "&fq=idxtype:TAXON&pageSize=50"
  )
  raw <- tryCatch(jsonlite::fromJSON(u, simplifyVector = FALSE, simplifyDataFrame = FALSE), error = function(e) e)
  if (inherits(raw, "error")) stop("NBN species-ws search failed: ", conditionMessage(raw))
  res <- raw$searchResults$results
  if (is.null(res) || length(res) == 0) return(data.frame())
  # flatten minimally
  out <- lapply(res, function(x) {
    data.frame(
      guid = x$guid %||% NA_character_,
      name = x$name %||% NA_character_,
      scientificName = x$scientificName %||% NA_character_,
      rank = x$rank %||% NA_character_,
      kingdom = x$kingdom %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  out
}

nbn_pick_guid <- function(results_df, sp) {
  if (!is.data.frame(results_df) || nrow(results_df) == 0) return(NA_character_)
  # Prefer exact match on scientificName or name
  sp_l <- tolower(sp)
  results_df <- results_df %>%
    mutate(
      sci_l = tolower(scientificName),
      name_l = tolower(name),
      exact = (sci_l == sp_l) | (name_l == sp_l)
    ) %>%
    arrange(desc(exact))
  guid <- results_df$guid[1]
  if (is.na(guid) || !nzchar(guid)) return(NA_character_)
  guid
}

# ---- records-ws totals + download --------------------------------------------

nbn_records_ws_total_for_fq <- function(guid, fq_vec) {
  u <- paste0(
    "https://records-ws.nbnatlas.org/occurrences/search?",
    nbn_build_query(list(
      q = paste0("lsid:", guid),
      fq = fq_vec,
      pageSize = 0
    ))
  )
  raw <- tryCatch(jsonlite::fromJSON(u), error = function(e) e)
  if (inherits(raw, "error")) return(NA_integer_)
  suppressWarnings(as.integer(raw$totalRecords))
}

download_zip <- function(url, dest) {
  old_timeout <- getOption("timeout")
  options(timeout = max(as.integer(old_timeout), 600L))
  on.exit(options(timeout = old_timeout), add = TRUE)
  
  if (requireNamespace("curl", quietly = TRUE)) {
    tryCatch({ curl::curl_download(url, destfile = dest, quiet = TRUE, mode = "wb"); TRUE },
             error = function(e) FALSE)
  } else {
    rc <- tryCatch(utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE, method = "libcurl"),
                   error = function(e) 1L, warning = function(w) 1L)
    if (!identical(rc, 0L)) {
      rc <- tryCatch(utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE),
                     error = function(e) 1L, warning = function(w) 1L)
    }
    identical(rc, 0L)
  }
}

nbn_records_ws_download_chunk <- function(guid, fq_vec, work_root, slug) {
  dl_base <- "https://records-ws.nbnatlas.org/occurrences/index/download"
  params <- list(
    q = paste0("lsid:", guid),
    fq = fq_vec,
    email = nbn_email,
    reasonTypeId = download_reason_id,
    fileType = "csv",
    dwcHeaders = "true",
    qa = "none"
  )
  dl_url <- paste0(dl_base, "?", nbn_build_query(params))
  
  mk_tag <- function() {
    paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)))
  }
  
  last_err <- NA_character_
  last_http <- NA_integer_
  
  for (attempt in seq_len(dl_max_tries)) {
    tag <- mk_tag()
    zip_path <- file.path(work_root, paste0("nbn_topup_", slug, "_", tag, ".zip"))
    unzip_dir <- file.path(work_root, paste0("nbn_topup_", slug, "_", tag, "_unzipped"))
    
    ok <- download_zip(dl_url, zip_path)
    
    if (!file.exists(zip_path) || is.na(file.info(zip_path)$size)) {
      last_http <- http_head_status(dl_url)
      last_err <- paste0("zip_missing (http=", last_http %||% NA_integer_, ")")
    } else {
      sz <- file.info(zip_path)$size
      
      if (looks_like_html(zip_path)) {
        last_http <- http_head_status(dl_url)
        last_err <- paste0("not_a_zip_html (bytes=", sz, "; http=", last_http %||% NA_integer_, ")")
        if (isTRUE(dl_debug_keep_failed)) {
          save_failed_payload_preview(
            zip_path,
            file.path(work_root, paste0("FAILED_", slug, "_", tag, ".html"))
          )
        }
      } else if (sz < dl_min_zip_bytes) {
        last_http <- http_head_status(dl_url)
        last_err <- paste0("zip_too_small (bytes=", sz, "; http=", last_http %||% NA_integer_, ")")
      } else {
        dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
        
        unz_ok <- tryCatch({
          utils::unzip(zip_path, exdir = unzip_dir)
          TRUE
        }, error = function(e) FALSE)
        
        if (!isTRUE(unz_ok)) {
          last_http <- http_head_status(dl_url)
          last_err <- paste0("unzip_failed (http=", last_http %||% NA_integer_, ")")
        } else {
          csv_path <- file.path(unzip_dir, "data.csv")
          if (!file.exists(csv_path)) {
            csvs <- list.files(unzip_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
            if (length(csvs) == 0) {
              last_err <- "no_csv_in_zip"
              csv_path <- NA_character_
            } else {
              csv_path <- csvs[which.max(file.info(csvs)$size)]
            }
          }
          
          if (!is.na(csv_path) && file.exists(csv_path)) {
            n2 <- tryCatch(length(readLines(csv_path, n = 2, warn = FALSE)), error = function(e) 0L)
            if (n2 >= 2) {
              return(list(
                ok = TRUE,
                zip_path = zip_path,
                csv_path = csv_path,
                unzip_dir = unzip_dir,
                error = NA_character_,
                http = last_http %||% NA_integer_,
                attempt = attempt,
                url = dl_url
              ))
            } else {
              last_err <- "csv_no_rows"
            }
          } else if (is.na(csv_path)) {
            last_err <- last_err %||% "csv_missing_after_unzip"
          } else {
            last_err <- "csv_missing_after_unzip"
          }
        }
      }
    }
    
    if (attempt < dl_max_tries) {
      sleep_s <- dl_backoff_base_sec * attempt
      message("[1.6]     Retry ", attempt, "/", dl_max_tries - 1L,
              " after ", sleep_s, "s (", last_err, ")")
      Sys.sleep(sleep_s)
    }
  }
  
  list(
    ok = FALSE,
    zip_path = NA_character_,
    csv_path = NA_character_,
    unzip_dir = NA_character_,
    error = paste0("Download failed after ", dl_max_tries, " tries: ", last_err),
    http = last_http %||% NA_integer_,
    attempt = dl_max_tries,
    url = dl_url
  )
}

# ---- Standardise records-ws CSV to Stage 1 NBN schema -------------------------

standardise_nbn_records_ws_df <- function(df, species_name) {
  # Case-insensitive column picker (keeps original name)
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
  col_basis  <- pick_col(c("basisOfRecord", "basis_of_record"))
  col_rank   <- pick_col(c("taxonRank", "taxon_rank", "rank"))
  col_occst  <- pick_col(c("occurrenceStatus", "occurrence_status"))
  col_dk     <- pick_col(c("datasetKey", "dataset_key"))
  col_dn     <- pick_col(c("datasetName", "dataset_name"))
  col_pok    <- pick_col(c("publishingOrgKey", "publishing_org_key"))
  col_inst   <- pick_col(c("institutionCode", "institution_code"))
  col_coll   <- pick_col(c("collectionCode", "collection_code"))
  
  recordID <- if (!is.null(col_record)) as.character(df[[col_record]]) else NA_character_
  scientificName <- if (!is.null(col_sci)) as.character(df[[col_sci]]) else species_name
  eventDate <- if (!is.null(col_date)) as.character(df[[col_date]]) else NA_character_
  eventDate <- as.character(eventDate)  # ensure always character (even if readr parsed Date)
  
  year <- if (!is.null(col_year)) {
    suppressWarnings(as.integer(df[[col_year]]))
  } else {
    suppressWarnings(as.integer(substr(as.character(eventDate), 1, 4)))
  }
  
  lat <- if (!is.null(col_lat)) suppressWarnings(as.numeric(df[[col_lat]])) else NA_real_
  lon <- if (!is.null(col_lon)) suppressWarnings(as.numeric(df[[col_lon]])) else NA_real_
  
  license <- if (!is.null(col_lic)) as.character(df[[col_lic]]) else NA_character_
  cuim <- if (!is.null(col_cuim)) suppressWarnings(as.numeric(df[[col_cuim]])) else NA_real_
  cp   <- if (!is.null(col_cp)) as.character(df[[col_cp]]) else NA_character_
  iv   <- if (!is.null(col_iv)) as.character(df[[col_iv]]) else NA_character_
  idby <- if (!is.null(col_idby)) as.character(df[[col_idby]]) else NA_character_
  
  basis <- if (!is.null(col_basis)) as.character(df[[col_basis]]) else NA_character_
  rank  <- if (!is.null(col_rank)) as.character(df[[col_rank]]) else NA_character_
  occst <- if (!is.null(col_occst)) as.character(df[[col_occst]]) else NA_character_
  dk    <- if (!is.null(col_dk)) as.character(df[[col_dk]]) else NA_character_
  dn    <- if (!is.null(col_dn)) as.character(df[[col_dn]]) else NA_character_
  pok   <- if (!is.null(col_pok)) as.character(df[[col_pok]]) else NA_character_
  inst  <- if (!is.null(col_inst)) as.character(df[[col_inst]]) else NA_character_
  coll  <- if (!is.null(col_coll)) as.character(df[[col_coll]]) else NA_character_
  
  out <- tibble(
    source = "NBN",
    species = species_name,
    recordID = recordID,
    lon = lon,
    lat = lat,
    date = eventDate,
    year = year,
    licence_raw = license,
    licence = license,                 # keep raw; Stage 3 can normalise further
    licence_expected = NA,             # Stage 1 had this; keep placeholder for schema consistency
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
    collectionCode = coll
  )
  
  out
}

# ---- Chunk planning: recursively split year ranges if near cap ----------------

near_cap <- function(total) {
  isTRUE(!is.na(total) && total >= (nbn_cap_n * nbn_near_cap_prop))
}

plan_year_ranges <- function(guid, y1, y2, depth = 0L) {
  if (depth > max_split_depth) {
    return(list(list(y1 = y1, y2 = y2, total = NA_integer_, status = "max_depth")))
  }
  
  fq_vec <- c(
    '-occurrence_status:"absent"',
    paste0("year:[", y1, " TO ", y2, "]")
  )
  
  tot <- nbn_records_ws_total_for_fq(guid, fq_vec)
  
  # If totals call fails, keep the chunk but mark unknown
  if (is.na(tot)) {
    return(list(list(y1 = y1, y2 = y2, total = NA_integer_, status = "total_unknown")))
  }
  
  # Empty chunk
  if (tot == 0L) {
    return(list(list(y1 = y1, y2 = y2, total = 0L, status = "empty")))
  }
  
  # Safe chunk
  if (!near_cap(tot)) {
    return(list(list(y1 = y1, y2 = y2, total = tot, status = "ok")))
  }
  
  # Near-cap
  if (y1 == y2) {
    # cannot split further by year in this version
    return(list(list(y1 = y1, y2 = y2, total = tot, status = "too_big_single_year")))
  }
  
  mid <- as.integer(floor((y1 + y2) / 2))
  left  <- plan_year_ranges(guid, y1, mid, depth + 1L)
  right <- plan_year_ranges(guid, mid + 1L, y2, depth + 1L)
  c(left, right)
}

# ---- Per-species checkpoint ---------------------------------------------------

read_state <- function(slug) {
  f <- file.path(nbn_topup_ckpt_dir, paste0("nbn_topup_state_", slug, ".rds"))
  if (!file.exists(f)) return(list(done = character(0), last_error = NA_character_, updated = NA_character_))
  tmp <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(tmp) || !is.list(tmp)) return(list(done = character(0), last_error = "state_read_error", updated = NA_character_))
  tmp
}

write_state <- function(slug, state) {
  f <- file.path(nbn_topup_ckpt_dir, paste0("nbn_topup_state_", slug, ".rds"))
  state$updated <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  saveRDS(state, f)
}

range_key <- function(y1, y2) paste0(y1, "_", y2)

# ---- Main loop ----------------------------------------------------------------

work_root <- file.path(nbn_topup_ckpt_dir, paste0("nbn_topup_work_", group_dir))
dir.create(work_root, recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()

for (i in seq_len(nrow(worklist))) {
  sp   <- worklist$species[i]
  slug <- worklist$slug[i] %||% slugify_species(sp)
  
  message("\n[1.6] (", i, "/", nrow(worklist), ") ", sp, " [", slug, "]")
  
  # Resolve GUID
  guid <- tryCatch({
    res <- nbn_species_ws_search(sp)
    nbn_pick_guid(res, sp)
  }, error = function(e) NA_character_)
  
  if (is.na(guid) || !nzchar(guid)) {
    message("[1.6]   GUID resolve FAILED; skipping.")
    summary_rows[[length(summary_rows) + 1L]] <- tibble(
      species = sp, slug = slug, guid = NA_character_,
      chunks_planned = NA_integer_, chunks_done = 0L,
      new_rows = 0L, final_rows = NA_integer_,
      status = "guid_failed"
    )
    next
  }
  
  # Plan year-range chunks
  chunks <- plan_year_ranges(guid, start_year, end_year)
  chunks_df <- bind_rows(lapply(chunks, as_tibble)) %>%
    mutate(key = range_key(y1, y2))
  
  # Remove empty chunks from plan (no need to download)
  chunks_df <- chunks_df %>% filter(!(status == "empty" & total == 0L))
  
  message("[1.6]   Planned chunks: ", nrow(chunks_df))
  if (nrow(chunks_df) > 0) {
    # show a short preview
    print(head(chunks_df %>% select(y1, y2, total, status), 10))
  }
  
  # Load state (what’s already done)
  state <- read_state(slug)
  done_keys <- unique(state$done %||% character(0))
  
  # Only attempt chunks not done
  todo <- chunks_df %>% filter(!key %in% done_keys)
  
  if (nrow(todo) == 0) {
    message("[1.6]   Nothing new to do (all chunks already done).")
    # still report file rows
    out_file <- file.path(nbn_run_dir, paste0("nbn_", slug, "_clean.csv"))
    final_rows <- if (file.exists(out_file)) {
      tryCatch(nrow(readr::read_csv(out_file, show_col_types = FALSE)), error = function(e) NA_integer_)
    } else NA_integer_
    
    summary_rows[[length(summary_rows) + 1L]] <- tibble(
      species = sp, slug = slug, guid = guid,
      chunks_planned = nrow(chunks_df), chunks_done = nrow(chunks_df),
      new_rows = 0L, final_rows = final_rows,
      status = "already_complete"
    )
    next
  }
  
  if (isTRUE(dry_run)) {
    message("[1.6]   dry_run=TRUE; skipping downloads.")
    summary_rows[[length(summary_rows) + 1L]] <- tibble(
      species = sp, slug = slug, guid = guid,
      chunks_planned = nrow(chunks_df), chunks_done = length(done_keys),
      new_rows = 0L, final_rows = NA_integer_,
      status = "dry_run"
    )
    next
  }
  
  # Load existing output (if present)
  out_file <- file.path(nbn_run_dir, paste0("nbn_", slug, "_clean.csv"))
  existing <- tryCatch(
    readr::read_csv(out_file, show_col_types = FALSE, col_types = cols(date = col_character())),
    error = function(e) NULL
  )
  
  new_rows_total <- 0L
  chunks_done_now <- 0L
  errors_now <- character(0)
  
  # Process todo chunks
  for (j in seq_len(nrow(todo))) {
    y1 <- todo$y1[j]
    y2 <- todo$y2[j]
    k  <- todo$key[j]
    st <- todo$status[j]
    tot <- todo$total[j]
    
    message("[1.6]   Download chunk ", j, "/", nrow(todo), ": years ", y1, "-", y2,
            " (status=", st, ", total=", tot %||% NA_integer_, ")")
    
    fq_vec <- c(
      '-occurrence_status:"absent"',
      paste0("year:[", y1, " TO ", y2, "]")
    )
    
    dl <- nbn_records_ws_download_chunk(guid, fq_vec, work_root = work_root, slug = slug)
    
    if (!isTRUE(dl$ok)) {
      msg_http <- if (!is.null(dl$http) && !is.na(dl$http)) paste0(" (http=", dl$http, ")") else ""
      message("[1.6]     FAILED: ", dl$error, msg_http)
      errors_now <- c(errors_now, paste0(k, ":", dl$error))
      state$last_error <- dl$error
      write_state(slug, state)
      next
    }
    
    # Read downloaded CSV and standardise
    raw_df <- tryCatch(readr::read_csv(dl$csv_path, show_col_types = FALSE), error = function(e) e)
    if (inherits(raw_df, "error")) {
      message("[1.6]     FAILED reading CSV: ", conditionMessage(raw_df))
      errors_now <- c(errors_now, paste0(k, ":read_csv"))
      state$last_error <- "read_csv"
      write_state(slug, state)
      next
    }
    
    std <- standardise_nbn_records_ws_df(raw_df, species_name = sp)
    
    if (!("recordID" %in% names(std))) {
      message("[1.6]     FAILED: recordID missing after standardisation.")
      errors_now <- c(errors_now, paste0(k, ":no_recordID"))
      state$last_error <- "no_recordID"
      write_state(slug, state)
      next
    }
    
    # Append into accumulator (force stable types to avoid bind_rows Date/character clashes)
    if (is.null(existing)) {
      existing <- std
    } else {
      if ("date" %in% names(existing)) existing$date <- as.character(existing$date)
      if ("date" %in% names(std))      std$date      <- as.character(std$date)
      
      if ("year" %in% names(existing)) existing$year <- suppressWarnings(as.integer(existing$year))
      if ("year" %in% names(std))      std$year      <- suppressWarnings(as.integer(std$year))
      
      existing <- bind_rows(existing, std)
    }
    
    new_rows_total <- new_rows_total + nrow(std)
    chunks_done_now <- chunks_done_now + 1L
    
    # Mark done in checkpoint
    state$done <- unique(c(done_keys, state$done %||% character(0), k))
    state$last_error <- NA_character_
    write_state(slug, state)
  }
  
  # If nothing appended, skip write
  if (is.null(existing)) {
    summary_rows[[length(summary_rows) + 1L]] <- tibble(
      species = sp, slug = slug, guid = guid,
      chunks_planned = nrow(chunks_df), chunks_done = length(unique(state$done %||% character(0))),
      new_rows = 0L, final_rows = NA_integer_,
      status = if (length(errors_now)) "failed_no_data" else "no_data"
    )
    next
  }
  
  # De-duplicate by recordID (keep first)
  before_n <- nrow(existing)
  existing <- existing %>%
    mutate(recordID = as.character(recordID)) %>%
    filter(!is.na(recordID) & nzchar(recordID)) %>%
    distinct(recordID, .keep_all = TRUE)
  after_n <- nrow(existing)
  
  # Write updated file
  readr::write_csv(existing, out_file, na = "")
  message("[1.6]   Wrote: ", normalizePath(out_file, winslash = "/", mustWork = FALSE),
          " (rows ", before_n, " -> ", after_n, " after dedupe)")
  
  # Summarise species status
  any_too_big <- any(chunks_df$status == "too_big_single_year", na.rm = TRUE)
  status <- if (length(errors_now) > 0) {
    if (any_too_big) "partial_with_too_big_year" else "partial_with_errors"
  } else {
    if (any_too_big) "complete_but_too_big_year" else "complete"
  }
  
  summary_rows[[length(summary_rows) + 1L]] <- tibble(
    species = sp, slug = slug, guid = guid,
    chunks_planned = nrow(chunks_df),
    chunks_done = length(unique(state$done %||% character(0))),
    new_rows = new_rows_total,
    final_rows = after_n,
    status = status
  )
}

# ---- Write one summary CSV ----------------------------------------------------

summary_df <- bind_rows(summary_rows)

dir.create(audits_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(summary_df, out_summary_csv, na = "")

message("\n[1.6] Wrote summary: ", normalizePath(out_summary_csv, winslash = "/", mustWork = FALSE))
message("[1.6] Done.")
