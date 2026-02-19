# stage_01_audit_raw_pull.R
#
# InfluentialSpecies — Stage 1 audit (read-only)
#
# Goal (what this is for)
#   Given a single raw pull run folder (group_dir), quickly answer:
#     - Did we produce the expected GBIF + NBN outputs for every species in the meta list?
#     - Are those outputs "merge-ready" for Stage 2 (i.e., have required columns)?
#     - Which species should Stage 1.6 (NBN top-up) attempt to repair / extend?
#
# Design principles
#   - One primary output file: a single "manifest" CSV (human-readable *and* machine-consumable for Stage 1.5).
#   - The first columns carry the headline status; details are kept to the right.
#   - Console prints a compact scoreboard (X/N ok, missing, empty, etc.).
#   - Extra columns in CSVs are allowed (we only enforce required columns).
#   - Read-only: this script never modifies raw pull outputs or checkpoints.
#
# Notes about storage
#   This script uses base R file listing/reading (file.exists/list.files/readr) so it works for:
#     - local filesystem
#     - mounted network drives
#     - many "cloud VM with mounted bucket" setups
#   If your raw outputs are in object storage without a filesystem mount (pure GCS API),
#   you’ll need a small adapter to list/read objects; the status logic below stays the same.
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

# The pull run folder you want to audit (under data/raw/gbif/ and data/raw/nbn/).
group_dir <- "home_run_true_list"

# Repo structure (edit these if you move folders in future).
repo_root <- getwd()
meta_dir  <- file.path(repo_root, "data", "_meta")
audits_dir <- file.path(meta_dir, "audits")

gbif_run_dir <- file.path(repo_root, "data", "raw", "gbif", group_dir)
nbn_run_dir  <- file.path(repo_root, "data", "raw", "nbn",  group_dir)

# Meta list (canonical): headerless binomial CSV (one species per line).
# This is the authoritative list used to name/slug Stage 0 outputs.
meta_species_file <- file.path(meta_dir, "species_list_binomial.csv")

# Kept for compatibility if you swap back to Excel later (ignored for the headerless CSV).
meta_sheet <- NULL
meta_species_col <- "species"

# IMPORTANT:
#   Our current Excel has blank header cells that readxl names as ...1, ...2, ...3, etc.
#   In your sheet, the authoritative binomial column is currently called "...4".
#   (Best long-term fix is to rename that header in Excel to e.g. "ScientificName".)
meta_species_col <- "Species"

# NBN "cap risk" probing (online). Turn off if you want a fast, offline-only audit.
do_nbn_online_probe <- TRUE

# NBN cap heuristics: used for flagging only (and for Stage 1.6 recommendations).
nbn_cap_n <- 500000L
nbn_near_cap_prop <- 0.98        # cap-risk if totalRecords >= cap*0.98
nbn_topup_floor_prop <- 0.90     # recommend top-up if totalRecords >= cap*0.90 (policy; tweak as desired)

# CSV output naming
timestamp_tag <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
out_manifest_csv <- file.path(audits_dir, paste0("audit_raw_pull_manifest_", group_dir, "_", timestamp_tag, ".csv"))

# ---- Required schemas (Stage 2 merge expectations) ---------------------------
# These are the "must have" columns for Stage 2 merge. Extra columns are OK.
# If you change the Stage 1 pull transmute outputs, update these lists to match.

required_cols_gbif <- c(
  "source","species","gbifID","occurrenceID","lon","lat","date","year","country",
  "licence_raw","licence","licence_expected",
  "coordinateUncertaintyInMeters","identificationVerificationStatus","issues","identifiedBy","dateIdentified",
  "basisOfRecord","taxonRank","occurrenceStatus","datasetKey","datasetName","publishingOrgKey","institutionCode","collectionCode"
)

required_cols_nbn <- c(
  "source","species","recordID","lon","lat","date","year",
  "licence_raw","licence","licence_expected",
  "coordinateUncertaintyInMeters","coordinatePrecision","identificationVerificationStatus","identifiedBy",
  "basisOfRecord","taxonRank","occurrenceStatus","datasetKey","datasetName","publishingOrgKey","institutionCode","collectionCode"
)

# ---- Helpers -----------------------------------------------------------------

stop_if_missing <- function(path, hint) {
  if (!dir.exists(path) && !file.exists(path)) {
    stop(
      paste0(
        "Missing path: ", path, "\n",
        hint
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

fmt_path <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

slugify_species <- function(species_name) {
  slug <- str_replace_all(tolower(species_name), "[^a-z0-9]+", "_")
  slug <- str_replace_all(slug, "^_+|_+$", "")
  slug
}

# Read only header (column names). Returns character vector or error marker.
read_csv_header <- function(path) {
  tryCatch(
    names(readr::read_csv(path, n_max = 0, show_col_types = FALSE, progress = FALSE)),
    error = function(e) structure(NA_character_, error = conditionMessage(e))
  )
}

# Cheap "header-only?" check without reading whole file:
# - reads up to 2 lines; if there's no 2nd line => empty_file (header-only).
is_header_only_csv <- function(path) {
  if (!file.exists(path)) return(NA)
  ln <- tryCatch(readLines(path, n = 2, warn = FALSE), error = function(e) character(0))
  length(ln) < 2
}

# Pull slug from a filename like gbif_<slug>_clean.csv / nbn_<slug>_clean.csv
slug_from_filename <- function(fname, prefix) {
  m <- str_match(fname, paste0("^", prefix, "_(.+)_clean\\.csv$"))
  if (is.na(m[1, 2])) return(NA_character_)
  m[1, 2]
}

# List all per-species output CSVs under a run directory (recursively).
# Returns a tibble with slug + path.
index_run_outputs <- function(run_dir, prefix) {
  if (!dir.exists(run_dir)) return(tibble(slug = character(0), path = character(0)))
  files <- list.files(run_dir, pattern = paste0("^", prefix, "_.+_clean\\.csv$"), recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(tibble(slug = character(0), path = character(0)))
  
  tibble(
    slug = vapply(basename(files), slug_from_filename, character(1), prefix = prefix),
    path = files
  ) %>%
    filter(!is.na(slug), nzchar(slug)) %>%
    group_by(slug) %>%
    # If multiple matches exist (shouldn't), keep the shortest path (closest to run root).
    arrange(nchar(path), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()
}

# Detect a meta list .xlsx if the default isn't present
find_meta_xlsx <- function(meta_dir) {
  xlsx <- list.files(meta_dir, pattern = "\\.xlsx$", full.names = TRUE)
  if (length(xlsx) == 0) return(NA_character_)
  xlsx[1]
}

# Read the meta species list from Excel.
#
# Important behaviour:
#   - This audit should *not* silently drop species rows.
#   - We only drop truly empty rows or rows we cannot parse into a binomial at all.
#
# Normalisation:
#   - If the cell contains "Common name (Genus species)" we extract the binomial in parentheses.
#   - Otherwise we take the first two tokens that look like a binomial from the start of the string.
read_meta_species <- function(xlsx_path, sheet = NULL, species_col = NULL) {
  
  if (!file.exists(xlsx_path)) {
    stop("Meta species file not found at: ", xlsx_path, call. = FALSE)
  }
  
  ext <- tolower(tools::file_ext(xlsx_path))
  
  # ---- Canonical binomial CSV (headerless, single column) ---------------------
  if (ext == "csv") {
    
    df <- readr::read_csv(
      xlsx_path,
      col_names = "species",
      show_col_types = FALSE,
      trim_ws = TRUE,
      progress = FALSE
    )
    
    sp <- df$species %>%
      as.character() %>%
      stringr::str_trim()
    
    sp <- sp[!is.na(sp) & nzchar(sp)]
    sp <- unique(sp)
    
    message("Meta species source used: ", basename(xlsx_path), " (headerless CSV)")
    message("Meta species count (unique): ", length(sp))
    
    return(sp)
  }
  
  # ---- Excel meta sheet (fallback / legacy) ----------------------------------
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to read the meta Excel file. Please install it.", call. = FALSE)
  }
  
  df <- readxl::read_excel(xlsx_path, sheet = sheet %||% 1)
  if (!is.data.frame(df) || nrow(df) == 0) stop("Meta Excel read returned no rows.", call. = FALSE)
  
  nm <- names(df)
  
  if (is.null(species_col)) {
    hit <- nm[str_detect(tolower(nm), "species")]
    species_col <- if (length(hit) == 0) nm[1] else hit[1]
  }
  
  if (!species_col %in% nm) stop("meta_species_col not found in meta sheet: ", species_col, call. = FALSE)
  
  sp_raw <- df[[species_col]] %>%
    as.character() %>%
    stringr::str_trim()
  
  sp_raw <- sp_raw[!is.na(sp_raw) & nzchar(sp_raw)]
  
  extract_binomial <- function(x) {
    m <- stringr::str_match(x, "\\(([A-Z][a-z-]+\\s+[a-z-]+)\\)")
    if (!is.na(m[, 2])) return(m[, 2])
    
    m2 <- stringr::str_match(x, "^([A-Z][a-z-]+)\\s+([a-z-]+)")
    if (!is.na(m2[, 1])) return(paste(m2[, 2], m2[, 3]))
    
    NA_character_
  }
  
  sp <- vapply(sp_raw, extract_binomial, character(1))
  
  bad <- is.na(sp)
  if (any(bad)) {
    warning(
      "Could not extract a binomial for ",
      sum(bad),
      " meta rows; they will be ignored. Examples: ",
      paste(utils::head(sp_raw[bad], 5), collapse = " | ")
    )
  }
  
  sp <- unique(sp[!bad])
  
  message("Meta species column used: ", species_col)
  message("Meta species count (unique): ", length(sp))
  
  sp
}

# ---- NBN web services (optional enrichment) ----------------------------------

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

nbn_pick_guid <- function(results_df, sp) {
  if (!is.data.frame(results_df) || nrow(results_df) == 0) return(NA_character_)
  
  res2 <- results_df %>%
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
  
  # pageSize=1 forces the service to compute totalRecords.
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

# ---- Status logic -------------------------------------------------------------

status_from_file <- function(path, required_cols) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(list(
      status = "missing_file",
      header_only = NA,
      missing_required_cols = NA_character_,
      extra_cols = NA_character_,
      read_error = NA_character_
    ))
  }
  
  header_only <- is_header_only_csv(path)
  
  cols <- read_csv_header(path)
  read_error <- attr(cols, "error", exact = TRUE)
  if (length(cols) == 1 && is.na(cols)) {
    return(list(
      status = "unreadable",
      header_only = header_only,
      missing_required_cols = NA_character_,
      extra_cols = NA_character_,
      read_error = as.character(read_error %||% "unknown error")
    ))
  }
  
  missing_req <- setdiff(required_cols, cols)
  extra <- setdiff(cols, required_cols)
  
  if (length(missing_req) > 0) {
    return(list(
      status = "missing_required_columns",
      header_only = header_only,
      missing_required_cols = paste(missing_req, collapse = ";"),
      extra_cols = if (length(extra) > 0) paste(extra, collapse = ";") else "",
      read_error = NA_character_
    ))
  }
  
  if (isTRUE(header_only)) {
    return(list(
      status = "empty_file",
      header_only = TRUE,
      missing_required_cols = "",
      extra_cols = if (length(extra) > 0) paste(extra, collapse = ";") else "",
      read_error = NA_character_
    ))
  }
  
  list(
    status = "ok",
    header_only = FALSE,
    missing_required_cols = "",
    extra_cols = if (length(extra) > 0) paste(extra, collapse = ";") else "",
    read_error = NA_character_
  )
}

# Stage 2 readiness is intentionally strict:
#   - GBIF must be ok (non-empty + required columns present).
#   - NBN must be ok OR empty_file (header-only but schema present).
# If you want to allow missing NBN for taxa with zero UK records, prefer making Stage 1 always
# write an empty schema-correct NBN CSV; then this logic remains stable.
is_stage2_ready <- function(gbif_status, nbn_status) {
  identical(gbif_status, "ok") && (nbn_status %in% c("ok", "empty_file"))
}

# "Top-up recommended" is a policy decision for Stage 1.6.
# Default rule:
#   - If NBN is merge-ready but NBN totalRecords suggests it is near/at the export cap, recommend top-up.
#   - If NBN is missing/unreadable/schema-missing, also recommend top-up (because Stage 2 will fail).
recommend_topup <- function(nbn_status, nbn_total_records) {
  if (nbn_status %in% c("missing_file", "unreadable", "missing_required_columns")) return(TRUE)
  if (is.na(nbn_total_records)) return(FALSE)  # no online evidence; don't auto-recommend
  nbn_total_records >= as.integer(nbn_cap_n * nbn_topup_floor_prop)
}

# ---- Run ---------------------------------------------------------------------

# Basic safety checks
stop_if_missing(meta_dir,  "Expected repo layout includes data/_meta/. Set repo_root correctly.")
dir.create(audits_dir, recursive = TRUE, showWarnings = FALSE)

# Find meta xlsx if the default path doesn't exist
if (!file.exists(meta_species_file)) {
  auto_xlsx <- find_meta_xlsx(meta_dir)
  if (!is.na(auto_xlsx)) {
    meta_species_file <- auto_xlsx
  }
}

# Read meta species list
meta_species <- read_meta_species(meta_species_file, sheet = meta_sheet, species_col = meta_species_col)

# Guardrail: if we unexpectedly read very few species, fail fast.
if (length(meta_species) < 10) {
  stop(
    "Meta species list looks too short (", length(meta_species), "). ",
    "Check meta_species_col and the contents of the Excel sheet."
  )
}

meta_tbl <- tibble(
  species = meta_species,
  slug = vapply(meta_species, slugify_species, character(1))
)

# Index run outputs (recursively) to tolerate subfolders under group_dir
gbif_idx <- index_run_outputs(gbif_run_dir, "gbif")
nbn_idx  <- index_run_outputs(nbn_run_dir,  "nbn")

message("Found outputs: GBIF ", nrow(gbif_idx), " | NBN ", nrow(nbn_idx))

# Find extras (outputs not in meta slug list)
meta_slugs <- meta_tbl$slug
gbif_extras <- gbif_idx %>% filter(!slug %in% meta_slugs)
nbn_extras  <- nbn_idx  %>% filter(!slug %in% meta_slugs)

# Join expected species to discovered files
manifest <- meta_tbl %>%
  left_join(gbif_idx %>% rename(gbif_path = path), by = "slug") %>%
  left_join(nbn_idx  %>% rename(nbn_path  = path), by = "slug") %>%
  mutate(
    gbif_path = ifelse(is.na(gbif_path), "", fmt_path(gbif_path)),
    nbn_path  = ifelse(is.na(nbn_path),  "", fmt_path(nbn_path))
  )

# Compute file/schema statuses
gbif_status_list <- lapply(manifest$gbif_path, status_from_file, required_cols = required_cols_gbif)
nbn_status_list  <- lapply(manifest$nbn_path,  status_from_file, required_cols = required_cols_nbn)

manifest <- manifest %>%
  mutate(
    gbif_status = vapply(gbif_status_list, `[[`, character(1), "status"),
    nbn_status  = vapply(nbn_status_list,  `[[`, character(1), "status"),
    
    gbif_header_only = vapply(gbif_status_list, function(x) as.logical(x$header_only), logical(1)),
    nbn_header_only  = vapply(nbn_status_list,  function(x) as.logical(x$header_only),  logical(1)),
    
    gbif_missing_required_cols = vapply(gbif_status_list, `[[`, character(1), "missing_required_cols"),
    nbn_missing_required_cols  = vapply(nbn_status_list,  `[[`, character(1), "missing_required_cols"),
    
    gbif_extra_cols = vapply(gbif_status_list, `[[`, character(1), "extra_cols"),
    nbn_extra_cols  = vapply(nbn_status_list,  `[[`, character(1), "extra_cols"),
    
    gbif_read_error = vapply(gbif_status_list, `[[`, character(1), "read_error"),
    nbn_read_error  = vapply(nbn_status_list,  `[[`, character(1), "read_error")
  ) %>%
  mutate(
    stage2_ready = mapply(is_stage2_ready, gbif_status, nbn_status) %>% as.logical()
  )

# Optional: enrich with NBN online totals and cap risk
if (isTRUE(do_nbn_online_probe)) {
  message("NBN online probe enabled (species-ws + records-ws). This can take a little while for ~100 species.")
  
  # Probe only where it makes sense (file exists and schema is merge-ready or empty).
  to_probe <- manifest %>%
    mutate(do_probe = nbn_status %in% c("ok", "empty_file")) %>%
    pull(do_probe)
  
  nbn_guid <- rep(NA_character_, nrow(manifest))
  nbn_total <- rep(NA_integer_, nrow(manifest))
  
  for (i in seq_len(nrow(manifest))) {
    if (!isTRUE(to_probe[i])) next
    sp <- manifest$species[i]
    
    res <- nbn_species_ws_search(sp)
    guid <- nbn_pick_guid(res, sp)
    nbn_guid[i] <- guid
    nbn_total[i] <- nbn_total_records(guid)
  }
  
  manifest <- manifest %>%
    mutate(
      nbn_guid = nbn_guid,
      nbn_total_records = nbn_total,
      nbn_cap_risk = !is.na(nbn_total_records) & (nbn_total_records >= as.integer(nbn_cap_n * nbn_near_cap_prop))
    )
} else {
  manifest <- manifest %>%
    mutate(
      nbn_guid = NA_character_,
      nbn_total_records = NA_integer_,
      nbn_cap_risk = NA
    )
}

# Decide which species Stage 1.6 should top up (policy-driven; easy to tweak)
manifest <- manifest %>%
  mutate(
    topup_recommended = mapply(recommend_topup, nbn_status, nbn_total_records) %>% as.logical(),
    topup_reason = case_when(
      nbn_status %in% c("missing_file") ~ "missing_nbn_file",
      nbn_status %in% c("unreadable") ~ "unreadable_nbn_file",
      nbn_status %in% c("missing_required_columns") ~ "nbn_missing_required_columns",
      !is.na(nbn_total_records) & nbn_total_records >= as.integer(nbn_cap_n * nbn_topup_floor_prop) ~ "nbn_near_cap_totalRecords",
      TRUE ~ ""
    ),
    needs_attention = (!stage2_ready) | isTRUE(topup_recommended) | isTRUE(nbn_cap_risk)
  )

# Reorder columns so the first few carry the main information.
# This table is intended to be:
#   - human-readable at a glance (first columns)
#   - directly usable as Stage 1.6 input (slug + topup_* fields)
manifest_out <- manifest %>%
  select(
    species,
    gbif_status,
    nbn_status,
    stage2_ready,
    topup_recommended,
    topup_reason,
    nbn_cap_risk,
    slug,
    gbif_path,
    nbn_path,
    gbif_missing_required_cols,
    nbn_missing_required_cols,
    gbif_read_error,
    nbn_read_error,
    gbif_extra_cols,
    nbn_extra_cols,
    nbn_total_records,
    nbn_guid,
    needs_attention
  )

# Write single manifest CSV
readr::write_csv(manifest_out, out_manifest_csv, na = "")
message("Wrote manifest: ", fmt_path(out_manifest_csv))

# ---- Console scoreboard -------------------------------------------------------

count_status <- function(x) as.list(table(factor(x, levels = c(
  "ok","empty_file","missing_file","missing_required_columns","unreadable"
))))

gbif_counts <- count_status(manifest_out$gbif_status)
nbn_counts  <- count_status(manifest_out$nbn_status)
N <- nrow(manifest_out)

stage2_ok_n <- sum(manifest_out$stage2_ready, na.rm = TRUE)
attn_n <- sum(manifest_out$needs_attention, na.rm = TRUE)
topup_n <- sum(manifest_out$topup_recommended, na.rm = TRUE)

gbif_extra_n <- nrow(gbif_extras)
nbn_extra_n  <- nrow(nbn_extras)

message("")
message("=== Stage 1.5 Audit Scoreboard (", group_dir, ") ===")
message("Meta species: ", N)
message("GBIF: ok ", gbif_counts$ok %||% 0, "/", N,
        " | missing ", gbif_counts$missing_file %||% 0,
        " | empty ", gbif_counts$empty_file %||% 0,
        " | missing_required_columns ", gbif_counts$missing_required_columns %||% 0,
        " | unreadable ", gbif_counts$unreadable %||% 0)

message("NBN : ok ", nbn_counts$ok %||% 0, "/", N,
        " | missing ", nbn_counts$missing_file %||% 0,
        " | empty ", nbn_counts$empty_file %||% 0,
        " | missing_required_columns ", nbn_counts$missing_required_columns %||% 0,
        " | unreadable ", nbn_counts$unreadable %||% 0)

message("Stage 2 merge-ready: ", stage2_ok_n, "/", N)

if (isTRUE(do_nbn_online_probe)) {
  cap_risk_n <- sum(manifest_out$nbn_cap_risk %in% TRUE, na.rm = TRUE)
  probed_n <- sum(manifest_out$nbn_status %in% c("ok","empty_file"))
  message("NBN cap-risk (online): ", cap_risk_n, " (probed ", probed_n, " species)")
} else {
  message("NBN cap-risk (online): [probe disabled]")
}

message("Top-up recommended (Stage 1.6 input): ", topup_n, " species")
message("Needs attention (any reason): ", attn_n, " species")
message("Extras on disk (not in meta): GBIF ", gbif_extra_n, " | NBN ", nbn_extra_n)

# If you want a quick peek at which species need attention, uncomment:
# print(manifest_out %>%
#         filter(needs_attention) %>%
#         select(species, gbif_status, nbn_status, stage2_ready, topup_recommended, topup_reason, nbn_cap_risk, slug))

# End of script
