# audit_raw_pull_v1_5.R
#
# InfluentialSpecies — Stage 1.5 audit (read-only)
#
# Purpose
#   Read-only QA over a raw pull run folder (GBIF + NBN) to:
#     1) confirm which species have outputs on disk (missing / extra vs the meta list)
#     2) verify output schemas are consistent (GBIF vs NBN; and across NBN pull routes)
#     3) flag NBN species that look at risk of truncation around the ~500k “download cap”
#        using NBN web services metadata (species-ws + records-ws totalRecords)
#
# Safety
#   - This script only READS raw CSVs + checkpoint folders.
#   - It does NOT overwrite raw outputs or checkpoints.
#   - It only writes NEW audit files under data/_meta/audits/.
#
# How to run
#   - Set your working directory to the repo root (folder containing data/), then:
#       source("R/audit_raw_pull_v1_5.R")
#   - Optionally edit the “User options” section below.
#
# Notes on NBN
#   - This audit uses NBN web services directly (species-ws + records-ws) for robust totals.
#   - It does not require galah_login().
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(jsonlite)
})

# ---- User options -------------------------------------------------------------

group_dir <- "home_run_2026-02-06"

# Meta list location (Excel). If this file is missing, the script will try to find an .xlsx in data/_meta/.
meta_dir <- file.path(getwd(), "data", "_meta")
meta_species_file <- file.path(meta_dir, "Influential Species Mapping List.xlsx")
meta_sheet <- NULL  # NULL = first sheet; or set e.g. "Mapping List"
meta_species_col <- NULL  # NULL = auto-detect; otherwise set to the column name containing species strings

# Expected NBN “cap” heuristic (used for flagging only)
nbn_cap_n <- 500000L
nbn_near_cap_prop <- 0.98  # “near cap” means >= 98% of cap
nbn_suspect_floor_prop <- 0.90  # do local row counting for totals >= 90% of cap

# NBN records-ws probe pacing (be gentle to the service)
nbn_pause_s <- 0.25

# Where to write audit outputs (new files only)
audit_out_dir <- file.path(meta_dir, "audits")

# If TRUE, attempt to count local NBN CSV rows for “suspect” species (reads the full CSV for those few species).
count_local_rows_for_suspects <- TRUE

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

fmt_path <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

cat_hr <- function(...) {
  cat("\n", paste(rep("=", 78), collapse = ""), "\n", sep = "")
  if (length(list(...)) > 0) cat(..., "\n", sep = "")
  invisible(TRUE)
}

stop_if_not_repo_root <- function(repo_root) {
  if (!dir.exists(file.path(repo_root, "data"))) {
    stop(
      "Working directory does not look like the repo root.\n",
      "Setwd to the InfluentialSpecies repo root (folder containing data/) and re-run.\n",
      "Current: ", fmt_path(repo_root)
    )
  }
  invisible(TRUE)
}

slugify_species <- function(species_name) {
  slug <- str_replace_all(tolower(species_name), "[^a-z0-9]+", "_")
  slug <- str_replace_all(slug, "^_+|_+$", "")
  slug
}

# Read only the header row (column names) from a CSV
read_csv_header_cols <- function(path) {
  if (!file.exists(path)) return(character())
  hdr <- tryCatch(readLines(path, n = 1, warn = FALSE, encoding = "UTF-8"), error = function(e) "")
  if (!nzchar(hdr)) return(character())
  strsplit(hdr, ",", fixed = TRUE)[[1]]
}

# Quick empty test: header-only CSVs are “empty”
csv_is_header_only <- function(path) {
  if (!file.exists(path)) return(NA)
  x <- tryCatch(readLines(path, n = 2, warn = FALSE, encoding = "UTF-8"), error = function(e) character())
  if (length(x) == 0) return(NA)
  length(x) == 1
}

# Auto-detect a likely meta .xlsx if the default file is missing
find_meta_xlsx <- function(meta_dir) {
  if (!dir.exists(meta_dir)) return(NA_character_)
  xlsx <- list.files(meta_dir, pattern = "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  if (length(xlsx) == 0) return(NA_character_)
  # Pick the most recently modified .xlsx
  fi <- file.info(xlsx)
  xlsx[order(fi$mtime, decreasing = TRUE)][1]
}

# Read expected species list (Latin binomials) from the mapping list sheet.
# The mapping list entries are typically like: "Brown bear (Ursus arctos)".
read_expected_species <- function(path, sheet = NULL, col = NULL) {
  if (!file.exists(path)) stop("Meta species file not found: ", fmt_path(path))
  
  # Prefer readxl if available; fall back to openxlsx.
  if (requireNamespace("readxl", quietly = TRUE)) {
    dat <- readxl::read_excel(path, sheet = sheet %||% 1)
    dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  } else if (requireNamespace("openxlsx", quietly = TRUE)) {
    sheet_use <- sheet %||% 1
    dat <- openxlsx::read.xlsx(path, sheet = sheet_use)
    dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  } else {
    stop("Need either readxl or openxlsx installed to read the meta .xlsx.")
  }
  
  if (nrow(dat) == 0) stop("Meta species file is readable but appears empty: ", fmt_path(path))
  
  nm <- names(dat)
  
  pick_col <- function(candidates) {
    hit <- intersect(candidates, nm)
    if (length(hit) == 0) NULL else hit[1]
  }
  
  col_use <- col
  if (is.null(col_use) || !nzchar(col_use) || !col_use %in% nm) {
    col_use <- pick_col(c("Species", "species", "ScientificName", "scientificName", "scientific_name"))
  }
  if (is.null(col_use)) stop("Could not auto-detect a species column in: ", fmt_path(path))
  
  raw <- dat[[col_use]]
  raw <- raw[!is.na(raw)]
  raw <- as.character(raw)
  raw <- str_trim(raw)
  raw <- raw[nzchar(raw)]
  
  # Extract the Latin name inside the last (...) group; otherwise keep if it already looks like a binomial.
  extract_latin <- function(x) {
    hits <- str_match_all(x, "\\(([^()]*)\\)")[[1]]
    if (nrow(hits) > 0) {
      latin <- str_trim(hits[nrow(hits), 2])
      return(latin)
    }
    # fallback: if it already looks like "Genus species"
    if (str_detect(x, "^[A-Z][a-z]+\\s+[a-z]+\\b")) {
      m <- str_match(x, "^([A-Z][a-z]+)\\s+([a-z]+)\\b")
      return(paste(m[2], m[3]))
    }
    x
  }
  
  sp <- vapply(raw, extract_latin, character(1))
  sp <- unique(sp)
  sp
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Detect whether species_subdir is being used for a given run root
detect_species_subdir <- function(root, prefix) {
  if (!dir.exists(root)) return(NA)
  flat <- list.files(root, pattern = paste0("^", prefix, "_.+_clean\\.csv$"), full.names = TRUE, recursive = FALSE)
  deep <- list.files(root, pattern = paste0("^", prefix, "_.+_clean\\.csv$"), full.names = TRUE, recursive = TRUE)
  if (length(flat) > 0) return(FALSE)
  if (length(deep) > 0) return(TRUE)
  NA
}

# Safe RDS read (read-only)
read_rds_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

# ------------------------------------------------------------------------------
# NBN web services: GUID resolution + totalRecords probe
# ------------------------------------------------------------------------------

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
  if (inherits(raw, "error")) return(data.frame())
  
  res <- raw$searchResults$results
  if (is.null(res) || length(res) == 0) return(data.frame())
  
  get1 <- function(x, nm, alt = NULL) {
    v <- x[[nm]]
    if ((is.null(v) || length(v) == 0) && !is.null(alt)) v <- x[[alt]]
    if (is.null(v) || length(v) == 0) return(NA_character_)
    as.character(v[[1]])
  }
  
  data.frame(
    scientificName  = vapply(res, get1, character(1), nm = "scientificName", alt = "name"),
    rank            = vapply(res, get1, character(1), nm = "rank"),
    taxonomicStatus = vapply(res, get1, character(1), nm = "taxonomicStatus"),
    guid            = vapply(res, get1, character(1), nm = "guid"),
    occurrenceCount = suppressWarnings(as.integer(vapply(res, get1, character(1), nm = "occurrenceCount"))),
    stringsAsFactors = FALSE
  )
}

nbn_pick_guid <- function(res, sp) {
  if (!is.data.frame(res) || nrow(res) == 0) return(NA_character_)
  
  res2 <- res %>%
    mutate(
      scientificName2 = tolower(as.character(scientificName)),
      rank2 = tolower(as.character(rank)),
      status2 = tolower(as.character(taxonomicStatus)),
      occ_n = suppressWarnings(as.integer(occurrenceCount))
    )
  
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

nbn_total_records <- function(guid) {
  if (is.na(guid) || !nzchar(guid)) return(NA_integer_)
  
  # pageSize=1 is used so the service definitely computes totalRecords.
  u <- paste0(
    "https://records-ws.nbnatlas.org/occurrences/search?",
    paste0(
      "q=", utils::URLencode(paste0("lsid:", guid), reserved = TRUE),
      "&fq=", utils::URLencode('-occurrence_status:"absent"', reserved = TRUE),
      "&pageSize=1&startIndex=0"
    )
  )
  
  raw <- tryCatch(jsonlite::fromJSON(u), error = function(e) e)
  if (inherits(raw, "error")) return(NA_integer_)
  
  suppressWarnings(as.integer(raw$totalRecords))
}

# Count lines efficiently for a small number of large CSVs (suspects only).
# Returns number of data rows (excludes header).
count_csv_rows_fast <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  
  if (requireNamespace("R.utils", quietly = TRUE)) {
    # countLines() counts all lines including header
    n_all <- tryCatch(R.utils::countLines(path), error = function(e) NA_integer_)
    if (is.na(n_all)) return(NA_integer_)
    return(max(0L, as.integer(n_all - 1L)))
  }
  
  # Fallback: chunked raw read counting '\n' (still reads the file, but avoids parsing CSV)
  con <- file(path, open = "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  
  buf_size <- 1024L * 1024L * 8L  # 8 MB
  n_newline <- 0L
  
  repeat {
    b <- readBin(con, what = "raw", n = buf_size)
    if (length(b) == 0) break
    n_newline <- n_newline + sum(b == as.raw(0x0A))
  }
  
  # Subtract 1 for header line if the file is non-empty
  max(0L, as.integer(n_newline - 1L))
}

# ------------------------------------------------------------------------------
# Main audit
# ------------------------------------------------------------------------------

repo_root <- getwd()
stop_if_not_repo_root(repo_root)

cat_hr("Stage 1.5 audit: starting")
cat("Repo root:  ", fmt_path(repo_root), "\n", sep = "")
cat("Run group:  ", group_dir, "\n", sep = "")
cat("Meta dir:   ", fmt_path(meta_dir), "\n", sep = "")

if (!file.exists(meta_species_file)) {
  alt <- find_meta_xlsx(meta_dir)
  if (!is.na(alt) && file.exists(alt)) {
    meta_species_file <- alt
  }
}
cat("Meta list:  ", fmt_path(meta_species_file), "\n", sep = "")

expected_species <- read_expected_species(meta_species_file, sheet = meta_sheet, col = meta_species_col)
expected_slugs <- vapply(expected_species, slugify_species, character(1))

# Run folders
gbif_run_root <- file.path(repo_root, "data", "raw", "gbif", group_dir)
nbn_run_root  <- file.path(repo_root, "data", "raw", "nbn",  group_dir)

cat_hr("Run folder checks")
cat("GBIF run root: ", fmt_path(gbif_run_root), " | exists=", dir.exists(gbif_run_root), "\n", sep = "")
cat("NBN  run root: ", fmt_path(nbn_run_root),  " | exists=", dir.exists(nbn_run_root),  "\n", sep = "")

if (!dir.exists(gbif_run_root) || !dir.exists(nbn_run_root)) {
  stop("Run folder(s) missing. Check group_dir and that you are setwd() to the repo root.")
}

gbif_species_subdir <- detect_species_subdir(gbif_run_root, "gbif")
nbn_species_subdir  <- detect_species_subdir(nbn_run_root,  "nbn")

cat_hr("Detected output layout")
cat("GBIF species_subdir guess: ", gbif_species_subdir, "\n", sep = "")
cat("NBN  species_subdir guess: ", nbn_species_subdir,  "\n", sep = "")

# Expected schemas (from the engine’s transmute outputs)
expected_cols_gbif <- c(
  "source","species","gbifID","occurrenceID","lon","lat","date","year","country",
  "licence_raw","licence","licence_expected",
  "coordinateUncertaintyInMeters","identificationVerificationStatus","issues","identifiedBy","dateIdentified",
  "basisOfRecord","taxonRank","occurrenceStatus","datasetKey","datasetName","publishingOrgKey","institutionCode","collectionCode"
)

expected_cols_nbn <- c(
  "source","species","recordID","lon","lat","date","year",
  "licence_raw","licence","licence_expected",
  "coordinateUncertaintyInMeters","coordinatePrecision","identificationVerificationStatus","identifiedBy",
  "basisOfRecord","taxonRank","occurrenceStatus","datasetKey","datasetName","publishingOrgKey","institutionCode","collectionCode"
)

# Build expected file paths for each species
build_outfile <- function(root, prefix, slug, species_subdir) {
  if (isTRUE(species_subdir)) {
    file.path(root, slug, paste0(prefix, "_", slug, "_clean.csv"))
  } else {
    file.path(root, paste0(prefix, "_", slug, "_clean.csv"))
  }
}

gbif_paths <- vapply(expected_slugs, function(slug) build_outfile(gbif_run_root, "gbif", slug, gbif_species_subdir), character(1))
nbn_paths  <- vapply(expected_slugs, function(slug) build_outfile(nbn_run_root,  "nbn",  slug, nbn_species_subdir),  character(1))

# Inventory the actual clean files present (for “extra files” detection)
gbif_present <- list.files(gbif_run_root, pattern = "^gbif_.+_clean\\.csv$", full.names = TRUE, recursive = isTRUE(gbif_species_subdir))
nbn_present  <- list.files(nbn_run_root,  pattern = "^nbn_.+_clean\\.csv$",  full.names = TRUE, recursive = isTRUE(nbn_species_subdir))

# Extract slug from filename
slug_from_file <- function(path, prefix) {
  b <- basename(path)
  m <- str_match(b, paste0("^", prefix, "_(.+)_clean\\.csv$"))
  if (is.na(m[1, 2])) NA_character_ else m[1, 2]
}

gbif_present_slugs <- vapply(gbif_present, slug_from_file, character(1), prefix = "gbif")
nbn_present_slugs  <- vapply(nbn_present,  slug_from_file, character(1), prefix = "nbn")

extras_gbif <- setdiff(gbif_present_slugs, expected_slugs)
extras_nbn  <- setdiff(nbn_present_slugs,  expected_slugs)

# Gather file stats + schema checks (header-only)
file_stats <- function(path, expected_cols) {
  ex <- file.exists(path)
  if (!ex) {
    return(list(
      exists = FALSE,
      bytes = NA_real_,
      mtime = NA_character_,
      header_only = NA,
      n_missing_cols = NA_integer_,
      n_extra_cols = NA_integer_,
      wrong_order = NA,
      missing_cols = NA_character_,
      extra_cols = NA_character_
    ))
  }
  fi <- file.info(path)
  cols <- read_csv_header_cols(path)
  
  miss <- setdiff(expected_cols, cols)
  extra <- setdiff(cols, expected_cols)
  wrong_order <- FALSE
  if (length(cols) == length(expected_cols) && all(sort(cols) == sort(expected_cols))) {
    wrong_order <- !identical(cols, expected_cols)
  }
  
  list(
    exists = TRUE,
    bytes = as.numeric(fi$size),
    mtime = as.character(fi$mtime),
    header_only = csv_is_header_only(path),
    n_missing_cols = length(miss),
    n_extra_cols = length(extra),
    wrong_order = wrong_order,
    missing_cols = if (length(miss) == 0) NA_character_ else paste(miss, collapse = ";"),
    extra_cols = if (length(extra) == 0) NA_character_ else paste(extra, collapse = ";")
  )
}

cat_hr("Scanning expected species (disk-only checks)")
audit <- tibble(
  species = expected_species,
  slug = expected_slugs,
  gbif_file = fmt_path(gbif_paths),
  nbn_file  = fmt_path(nbn_paths)
)

gbif_meta <- lapply(gbif_paths, file_stats, expected_cols = expected_cols_gbif)
nbn_meta  <- lapply(nbn_paths,  file_stats, expected_cols = expected_cols_nbn)

audit <- audit %>%
  mutate(
    gbif_exists = vapply(gbif_meta, `[[`, logical(1), "exists"),
    gbif_bytes  = vapply(gbif_meta, `[[`, numeric(1),  "bytes"),
    gbif_mtime  = vapply(gbif_meta, `[[`, character(1),"mtime"),
    gbif_header_only = vapply(gbif_meta, `[[`, logical(1),"header_only"),
    gbif_missing_cols_n = vapply(gbif_meta, `[[`, integer(1),"n_missing_cols"),
    gbif_extra_cols_n   = vapply(gbif_meta, `[[`, integer(1),"n_extra_cols"),
    gbif_wrong_order    = vapply(gbif_meta, `[[`, logical(1),"wrong_order"),
    gbif_missing_cols   = vapply(gbif_meta, `[[`, character(1),"missing_cols"),
    gbif_extra_cols     = vapply(gbif_meta, `[[`, character(1),"extra_cols"),
    
    nbn_exists = vapply(nbn_meta, `[[`, logical(1), "exists"),
    nbn_bytes  = vapply(nbn_meta, `[[`, numeric(1),  "bytes"),
    nbn_mtime  = vapply(nbn_meta, `[[`, character(1),"mtime"),
    nbn_header_only = vapply(nbn_meta, `[[`, logical(1),"header_only"),
    nbn_missing_cols_n = vapply(nbn_meta, `[[`, integer(1),"n_missing_cols"),
    nbn_extra_cols_n   = vapply(nbn_meta, `[[`, integer(1),"n_extra_cols"),
    nbn_wrong_order    = vapply(nbn_meta, `[[`, logical(1),"wrong_order"),
    nbn_missing_cols   = vapply(nbn_meta, `[[`, character(1),"missing_cols"),
    nbn_extra_cols     = vapply(nbn_meta, `[[`, character(1),"extra_cols")
  )

# Checkpoint folder sanity (read-only)
ckpt_root <- file.path(repo_root, "data", "_checkpoints")
gbif_ckpt_dir <- file.path(ckpt_root, "gbif")
nbn_ckpt_dir  <- file.path(ckpt_root, "nbn")

cat_hr("Checkpoint folder snapshot (read-only)")
cat("Checkpoint root: ", fmt_path(ckpt_root), "\n", sep = "")
cat("GBIF ckpt dir:   ", fmt_path(gbif_ckpt_dir), " | exists=", dir.exists(gbif_ckpt_dir), "\n", sep = "")
cat("NBN  ckpt dir:   ", fmt_path(nbn_ckpt_dir),  " | exists=", dir.exists(nbn_ckpt_dir),  "\n", sep = "")

n_gbif_ckpt <- if (dir.exists(gbif_ckpt_dir)) length(list.files(gbif_ckpt_dir, pattern = "^gbif_pull_checkpoint_.*\\.rds$", ignore.case = TRUE)) else 0L
n_nbn_state <- if (dir.exists(nbn_ckpt_dir))  length(list.files(nbn_ckpt_dir,  pattern = "^nbn_state_.*\\.rds$", ignore.case = TRUE)) else 0L
n_nbn_legacy <- if (dir.exists(nbn_ckpt_dir)) length(list.files(nbn_ckpt_dir, pattern = "^nbn_pull_checkpoint_.*\\.rds$", ignore.case = TRUE)) else 0L

cat("GBIF checkpoint files (gbif_pull_checkpoint_*): ", n_gbif_ckpt, "\n", sep = "")
cat("NBN state files (nbn_state_*):                 ", n_nbn_state, "\n", sep = "")
cat("NBN legacy files (nbn_pull_checkpoint_*):      ", n_nbn_legacy, "\n", sep = "")

# Attach GBIF checkpoint status when present
gbif_ckpt_paths <- vapply(expected_slugs, function(slug) file.path(gbif_ckpt_dir, paste0("gbif_pull_checkpoint_", slug, ".rds")), character(1))
gbif_ckpt_obj <- lapply(gbif_ckpt_paths, read_rds_safe)

get_ckpt_field <- function(x, field) {
  if (is.null(x) || !is.list(x) || is.null(x[[field]])) return(NA)
  x[[field]]
}

audit <- audit %>%
  mutate(
    gbif_ckpt_exists = file.exists(gbif_ckpt_paths),
    gbif_ckpt_complete = vapply(gbif_ckpt_obj, get_ckpt_field, logical(1), field = "complete"),
    gbif_ckpt_mode     = vapply(gbif_ckpt_obj, get_ckpt_field, character(1), field = "mode"),
    gbif_ckpt_expected = vapply(gbif_ckpt_obj, get_ckpt_field, integer(1), field = "total_expected"),
    gbif_ckpt_key      = vapply(gbif_ckpt_obj, get_ckpt_field, character(1), field = "download_key"),
    gbif_ckpt_status   = vapply(gbif_ckpt_obj, get_ckpt_field, character(1), field = "download_status"),
    gbif_ckpt_last_updated = vapply(gbif_ckpt_obj, get_ckpt_field, character(1), field = "last_updated")
  )

# NBN records-ws totals (service probe; read-only)
cat_hr("NBN records-ws totals (probe)")
cat("Probing totalRecords via species-ws GUID resolution + records-ws search.\n")
cat("Pause between calls: ", nbn_pause_s, "s\n", sep = "")

nbn_guid <- rep(NA_character_, length(expected_species))
nbn_total <- rep(NA_integer_, length(expected_species))
nbn_probe_note <- rep(NA_character_, length(expected_species))

for (i in seq_along(expected_species)) {
  sp <- expected_species[i]
  Sys.sleep(nbn_pause_s)
  
  res <- nbn_species_ws_search(sp)
  guid <- nbn_pick_guid(res, sp)
  nbn_guid[i] <- guid
  
  if (is.na(guid) || !nzchar(guid)) {
    nbn_probe_note[i] <- "no_exact_species_match_in_species_ws"
    next
  }
  
  Sys.sleep(nbn_pause_s)
  tot <- nbn_total_records(guid)
  nbn_total[i] <- tot
  
  if (is.na(tot)) {
    nbn_probe_note[i] <- "records_ws_totalRecords_failed"
  } else {
    nbn_probe_note[i] <- "ok"
  }
  
  if (i %% 10 == 0) cat("  probed ", i, " / ", length(expected_species), "\n", sep = "")
}

audit <- audit %>%
  mutate(
    nbn_guid = nbn_guid,
    nbn_totalRecords = nbn_total,
    nbn_probe_note = nbn_probe_note
  )

# Determine which species look “suspect” for local row counting (cap risk or inconsistencies)
audit <- audit %>%
  mutate(
    nbn_has_records_online = !is.na(nbn_totalRecords) & nbn_totalRecords > 0L,
    nbn_cap_risk_online = !is.na(nbn_totalRecords) & nbn_totalRecords >= nbn_cap_n,
    nbn_near_cap_online = !is.na(nbn_totalRecords) & nbn_totalRecords >= as.integer(floor(nbn_cap_n * nbn_near_cap_prop)),
    nbn_suspect_for_count = isTRUE(count_local_rows_for_suspects) & (
      (!is.na(nbn_totalRecords) & nbn_totalRecords >= as.integer(floor(nbn_cap_n * nbn_suspect_floor_prop))) |
        (nbn_has_records_online & (is.na(nbn_exists) | !nbn_exists | isTRUE(nbn_header_only))) |
        (nbn_cap_risk_online)
    )
  )

# Count local rows for suspect NBN species (read-only; reads file contents for those few)
nbn_local_rows <- rep(NA_integer_, nrow(audit))

if (isTRUE(count_local_rows_for_suspects)) {
  suspects <- which(audit$nbn_suspect_for_count)
  cat_hr("Counting local rows for suspect NBN CSVs")
  cat("Suspects to count: ", length(suspects), " / ", nrow(audit), "\n", sep = "")
  
  if (length(suspects) > 0) {
    for (k in seq_along(suspects)) {
      i <- suspects[k]
      path <- audit$nbn_file[i]
      if (!file.exists(path)) next
      cat("  [", k, "/", length(suspects), "] ", audit$species[i], " | counting rows...\n", sep = "")
      nbn_local_rows[i] <- count_csv_rows_fast(path)
    }
  }
}

audit <- audit %>%
  mutate(
    nbn_local_rows = nbn_local_rows,
    nbn_local_near_500k = !is.na(nbn_local_rows) & nbn_local_rows >= (nbn_cap_n - 2000L),
    nbn_local_vs_online_gap = ifelse(!is.na(nbn_local_rows) & !is.na(nbn_totalRecords),
                                     as.integer(nbn_totalRecords - nbn_local_rows),
                                     NA_integer_),
    nbn_possible_truncation = !is.na(nbn_totalRecords) & nbn_totalRecords > nbn_cap_n &
      !is.na(nbn_local_rows) & nbn_local_rows >= (nbn_cap_n - 2000L),
    nbn_possible_incomplete = nbn_has_records_online &
      (is.na(nbn_local_rows) | nbn_local_rows < pmax(0L, nbn_totalRecords - 1000L)) &
      (is.na(nbn_totalRecords) | nbn_totalRecords <= nbn_cap_n)  # “incomplete” without invoking the cap
  )

# Summary flags for downstream readiness (schema + presence)
audit <- audit %>%
  mutate(
    gbif_schema_ok = gbif_exists & !isTRUE(gbif_header_only) & gbif_missing_cols_n == 0L & gbif_extra_cols_n == 0L,
    nbn_schema_ok  = nbn_exists  & (!isTRUE(nbn_header_only) | (!is.na(nbn_totalRecords) & nbn_totalRecords == 0L)) &
      nbn_missing_cols_n == 0L & nbn_extra_cols_n == 0L,
    both_present = gbif_exists & nbn_exists,
    both_nonempty_or_valid = gbif_exists & !isTRUE(gbif_header_only) &
      (
        (!isTRUE(nbn_header_only) & nbn_exists) |
          (!is.na(nbn_totalRecords) & nbn_totalRecords == 0L)
      )
  )

# Detect missing/extra species relative to the meta list
missing_gbif <- audit %>% filter(!gbif_exists | isTRUE(gbif_header_only))
missing_nbn  <- audit %>% filter(!nbn_exists | (isTRUE(nbn_header_only) & (is.na(nbn_totalRecords) | nbn_totalRecords > 0L)))

# “Extra” files present in run folder that are not in the meta list
extras_tbl <- tibble(
  source = c(rep("GBIF", length(extras_gbif)), rep("NBN", length(extras_nbn))),
  slug = c(extras_gbif, extras_nbn)
)

# Run progress snapshot: newest file mtimes
latest_file_info <- function(files) {
  if (length(files) == 0) return(list(path = NA_character_, mtime = NA_character_))
  fi <- file.info(files)
  i <- which.max(fi$mtime)
  list(path = fmt_path(files[i]), mtime = as.character(fi$mtime[i]))
}

gbif_latest <- latest_file_info(gbif_present)
nbn_latest  <- latest_file_info(nbn_present)

# ------------------------------------------------------------------------------
# Write audit outputs (new files only)
# ------------------------------------------------------------------------------

dir.create(audit_out_dir, recursive = TRUE, showWarnings = FALSE)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

audit_file <- file.path(audit_out_dir, paste0("audit_raw_pull_", group_dir, "_", ts, ".csv"))
missing_file <- file.path(audit_out_dir, paste0("audit_missing_", group_dir, "_", ts, ".csv"))
cap_file <- file.path(audit_out_dir, paste0("audit_nbn_cap_risk_", group_dir, "_", ts, ".csv"))
extras_file <- file.path(audit_out_dir, paste0("audit_extras_", group_dir, "_", ts, ".csv"))
meta_sig_file <- file.path(audit_out_dir, paste0("audit_meta_signature_", group_dir, "_", ts, ".txt"))

write_csv(audit, audit_file)

missing_summary <- bind_rows(
  missing_gbif %>% transmute(source = "GBIF", species, slug, file = gbif_file),
  missing_nbn  %>% transmute(source = "NBN",  species, slug, file = nbn_file)
) %>% distinct()

write_csv(missing_summary, missing_file)

cap_risk <- audit %>%
  filter(nbn_near_cap_online | nbn_cap_risk_online | nbn_possible_truncation | nbn_possible_incomplete) %>%
  arrange(desc(nbn_totalRecords))

write_csv(cap_risk, cap_file)

write_csv(extras_tbl, extras_file)

meta_fi <- file.info(meta_species_file)
meta_sig <- c(
  paste0("meta_file=", fmt_path(meta_species_file)),
  paste0("meta_mtime=", as.character(meta_fi$mtime)),
  paste0("meta_bytes=", as.numeric(meta_fi$size)),
  paste0("meta_md5=", as.character(tools::md5sum(meta_species_file))),
  paste0("n_expected_species=", length(expected_species))
)
writeLines(meta_sig, meta_sig_file)

# ------------------------------------------------------------------------------
# Console summary
# ------------------------------------------------------------------------------

cat_hr("Audit summary")
cat("Expected species (from meta): ", length(expected_species), "\n", sep = "")
cat("GBIF clean files found in run: ", length(gbif_present), "\n", sep = "")
cat("NBN  clean files found in run: ", length(nbn_present),  "\n", sep = "")

cat("\nLatest GBIF file mtime: ", gbif_latest$mtime, "\n", sep = "")
cat("Latest GBIF file path:  ", gbif_latest$path,  "\n", sep = "")
cat("\nLatest NBN  file mtime: ", nbn_latest$mtime, "\n", sep = "")
cat("Latest NBN  file path:  ", nbn_latest$path,  "\n", sep = "")

cat_hr("Schema issues (header-only checks)")
cat("GBIF files with missing cols: ", sum(audit$gbif_missing_cols_n > 0, na.rm = TRUE), "\n", sep = "")
cat("GBIF files with extra cols:   ", sum(audit$gbif_extra_cols_n > 0,   na.rm = TRUE), "\n", sep = "")
cat("NBN  files with missing cols: ", sum(audit$nbn_missing_cols_n > 0,  na.rm = TRUE), "\n", sep = "")
cat("NBN  files with extra cols:   ", sum(audit$nbn_extra_cols_n > 0,    na.rm = TRUE), "\n", sep = "")

cat_hr("Missing/empty outputs (relative to meta list)")
cat("Missing/empty GBIF: ", nrow(missing_gbif), "\n", sep = "")
cat("Missing/empty NBN:  ", nrow(missing_nbn),  "\n", sep = "")

if (nrow(missing_summary) > 0) {
  cat("\nTop missing (first 20 rows):\n")
  print(missing_summary %>% slice_head(n = 20), n = 20)
}

cat_hr("NBN cap/truncation risk flags")
cat("Near-cap online totals (>= ", as.integer(floor(nbn_cap_n * nbn_near_cap_prop)), "): ",
    sum(audit$nbn_near_cap_online, na.rm = TRUE), "\n", sep = "")
cat("Online totals >= ", nbn_cap_n, ": ",
    sum(audit$nbn_cap_risk_online, na.rm = TRUE), "\n", sep = "")
cat("Possible truncation (online > cap AND local ~ cap): ",
    sum(audit$nbn_possible_truncation, na.rm = TRUE), "\n", sep = "")
cat("Possible incomplete (online > local but not invoking cap): ",
    sum(audit$nbn_possible_incomplete, na.rm = TRUE), "\n", sep = "")

if (nrow(cap_risk) > 0) {
  cat("\nTop cap/incomplete flags (first 20 rows):\n")
  print(
    cap_risk %>%
      transmute(
        species, slug,
        nbn_totalRecords,
        nbn_local_rows,
        nbn_possible_truncation,
        nbn_possible_incomplete,
        nbn_file
      ) %>%
      slice_head(n = 20),
    n = 20
  )
}

cat_hr("Extra files present (in run folder but not in meta list)")
cat("Extra GBIF slugs: ", length(extras_gbif), "\n", sep = "")
cat("Extra NBN  slugs: ", length(extras_nbn),  "\n", sep = "")

if (nrow(extras_tbl) > 0) {
  print(extras_tbl, n = nrow(extras_tbl))
}

cat_hr("Wrote audit outputs")
cat("Main audit CSV:      ", fmt_path(audit_file), "\n", sep = "")
cat("Missing summary CSV: ", fmt_path(missing_file), "\n", sep = "")
cat("Cap-risk summary CSV:", fmt_path(cap_file), "\n", sep = "")
cat("Extras summary CSV:  ", fmt_path(extras_file), "\n", sep = "")
cat("Meta signature TXT:  ", fmt_path(meta_sig_file), "\n", sep = "")
cat("\nDone.\n")
