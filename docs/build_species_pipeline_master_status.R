# docs/build_species_pipeline_master_status.R -----------------------------------
#
# InfluentialSpecies — Master species pipeline + uncertainty/obscuring status
#
# Purpose
#   Produce ONE sharing-friendly table that answers:
#     1) How far each species has progressed (downloaded → merged → filtered → gridded)
#     2) How many records we have at key steps (GBIF / NBN / merged / filtered / gridded)
#     3) How many records are lost to:
#          - cross-source de-duplication (GBIF vs NBN)
#          - policy filtering (overall)
#          - the uncertainty rule specifically (pre-filter, uncertainty-only)
#     4) For sensitive species, whether UK records show strong generalisation signals
#        (e.g., UK white-tailed eagle often reports ~7071 m uncertainty, which equals
#        5 km * sqrt(2) — the distance from the centre of a 10 km square to a corner).
#
# Key field
#   coordinateUncertaintyInMeters
#     A radius (meters) describing uncertainty around the reported coordinate.
#     Elevated values can reflect coarse coordinates, atlas squares, or intentional generalisation.
#
# Inputs
#   1) Canonical species list (defines the row set):
#        data/_meta/species_list_binomial.csv
#      Notes:
#        - This file has NO header.
#        - One binomial per line.
#
#   2) Runlogs + summary tables for pipeline counts (auto-detected / canonical-preferred):
#        data/processed/02_merged/_runlog_02_merged.csv (or 01_merged equivalent)
#        data/processed/04_filtered/_runlog_04_filtered.csv (or 03_filtered equivalent)
#        data/processed/**/_summary_grid.csv (latest found)
#
#   3) Optional sensitive-species list (for flagging + optional UK-only stats):
#        data/_meta/Combined-Sensitive-Species-List_06-25(.csv)
#
#   4) Occurrence outputs (for uncertainty audit; Stage 02 merged):
#        data/processed/02_merged/<slug>/occ_<slug>__merged.parquet
#
# Outputs
#   docs/derived/species_pipeline_master_status.csv
#   docs/derived/species_pipeline_master_status__README.txt
#
# Notes
#   - This script is intended for sharing (not diagnosis). It avoids per-stage status/note columns
#     and focuses on volumes + explainable losses.
#   - The master table row set is ALWAYS the canonical species list in data/_meta,
#     not whatever happens to exist in runlogs.
#   - "Uncertainty-rule drop" is computed by applying ONLY the uncertainty rule to Stage 02
#     merged records (ignoring other policy rules). This isolates obscuring/precision effects.
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# Uncertainty rule (must match the Stage 04 policy you are summarising)
default_max_uncertainty_m <- 1000L
uncertainty_missing_action <- "keep"   # "keep" or "drop"

# Sensitive list (used for flagging + UK-only stats for sensitive species)
# The file may exist with or without ".csv" — the script will find either.
sensitive_list_stem <- file.path("data", "_meta", "Combined-Sensitive-Species-List_06-25")

# UK-only stats are computed ONLY for sensitive species
compute_uk_stats_for_sensitive <- TRUE
uk_country_values <- c("United Kingdom", "UK", "GB")

# When runlogs are incomplete (e.g., partial reruns), backfill counts from files on disk.
backfill_counts_from_files <- TRUE

# Output paths (under docs/derived)
out_csv <- file.path("docs", "derived", "species_pipeline_master_status.csv")
write_readme <- TRUE

# ==============================================================================
# Repo root (robust when sourced from docs/ or scripts/)
# ==============================================================================

this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

find_repo_root <- function(start_dir) {
  markers <- c(".git", "R", "data", "InfluentialSpecies.Rproj", "DESCRIPTION")
  d <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:20) {
    if (any(file.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root <- if (!is.null(this_file) && nzchar(this_file)) {
  find_repo_root(dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE)))
} else {
  find_repo_root(getwd())
}

processed_root <- file.path(repo_root, "data", "processed")

# ==============================================================================
# Helpers
# ==============================================================================

yesno <- function(x) {
  x <- as.logical(x)
  ifelse(is.na(x), "", ifelse(x, "Yes", "No"))
}

slugify_species <- function(x) {
  s <- gsub("[^a-z0-9]+", "_", tolower(as.character(x)))
  gsub("^_+|_+$", "", s)
}

as_int <- function(x) suppressWarnings(as.integer(x))

pct1 <- function(num, den) {
  out <- rep(NA_real_, length(num))
  ok <- !is.na(num) & !is.na(den) & den > 0
  out[ok] <- round(100 * num[ok] / den[ok], 1)
  out
}

mode_numeric <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  tt <- sort(table(x), decreasing = TRUE)
  as.numeric(names(tt)[1])
}

safe_read_csv_dt <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  tryCatch(fread(file = path, fill = TRUE, na.strings = c("", "NA")), error = function(e) NULL)
}

read_parquet_cols <- function(path, cols) {
  tryCatch(arrow::read_parquet(path, col_select = cols), error = function(e) NULL)
}

norm_colnames <- function(nm) {
  nm <- tolower(as.character(nm))
  nm <- sub("^\ufeff", "", nm)               # strip BOM if present
  nm <- gsub("[^a-z0-9]+", "", nm)
  nm
}

pick_species_col <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  n0 <- names(dt)
  n1 <- norm_colnames(n0)
  
  # Common candidates (normalised)
  want <- c(
    "scientificname",
    "speciesbinomial",
    "binomial",
    "scientificnameaccepted",
    "taxon",
    "species",
    "taxonname",
    "acceptedname",
    "acceptedscientificname"
  )
  
  hit <- which(n1 %in% want)
  if (length(hit) > 0) return(n0[hit[1]])
  
  # Fall back to first column if nothing matches
  n0[1]
}

# Find a file either exactly (preferred) or by stem (with any extension)
find_by_stem <- function(stem_path_abs) {
  if (file.exists(stem_path_abs)) return(stem_path_abs)
  if (file.exists(paste0(stem_path_abs, ".csv"))) return(paste0(stem_path_abs, ".csv"))
  
  dir0 <- dirname(stem_path_abs)
  base0 <- basename(stem_path_abs)
  
  if (!dir.exists(dir0)) return(NA_character_)
  hits <- list.files(
    dir0,
    pattern = paste0("^", gsub("([\\.^$|()\\[\\]{}*+?\\\\-])", "\\\\\\1", base0), "(\\..*)?$"),
    full.names = TRUE,
    ignore.case = FALSE
  )
  hits <- hits[file.exists(hits)]
  if (length(hits) == 0) return(NA_character_)
  
  info <- file.info(hits)
  hits[which.max(info$mtime)]
}

find_latest_by_basename_regex <- function(root_dir, basename_regex_vec) {
  if (!dir.exists(root_dir)) return(NA_character_)
  all_files <- list.files(root_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  if (length(all_files) == 0) return(NA_character_)
  bn <- basename(all_files)
  
  keep <- rep(FALSE, length(all_files))
  for (rx in basename_regex_vec) keep <- keep | grepl(rx, bn, ignore.case = TRUE)
  
  hits <- all_files[keep]
  hits <- hits[file.exists(hits)]
  if (length(hits) == 0) return(NA_character_)
  
  info <- file.info(hits)
  hits[which.max(info$mtime)]
}

# Read one binomial per line from a headerless file
read_binomial_list_noheader <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(character())
  x <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- x[tolower(x) != "binomial"]
  unique(x)
}

# Normalise incoming text and extract a binomial from messy scientific names.
# This helps when lists include subspecies, aggregates, or short notes (e.g., "s. str.", "subsp. ...").
clean_text <- function(x) {
  x <- as.character(x)
  x <- sub("^\ufeff", "", x)                    # BOM
  x <- gsub("\u00A0", " ", x, fixed = TRUE)     # NBSP
  x <- gsub("\u200B|\u200C|\u200D|\uFEFF", "", x, perl = TRUE)  # zero-width
  x <- trimws(x)
  gsub("\\s+", " ", x)
}

extract_binomial <- function(x) {
  x <- clean_text(x)
  y <- gsub("[^A-Za-z\\-\\s]", " ", x)
  y <- gsub("\\s+", " ", trimws(y))
  if (!nzchar(y)) return(NA_character_)
  parts <- strsplit(y, " ", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) < 2) return(NA_character_)
  paste(parts[1], parts[2])
}

# ==============================================================================
# Canonical species list (defines the row set)
# ==============================================================================

species_list_path <- file.path(repo_root, "data", "_meta", "species_list_binomial.csv")
if (!file.exists(species_list_path)) stop("Can't find species list at: ", species_list_path)

species_binomials <- read_binomial_list_noheader(species_list_path)
if (length(species_binomials) == 0) stop("Species list read as 0 rows. Check: ", species_list_path)

base <- data.table(
  species = species_binomials,
  slug = slugify_species(species_binomials)
)

if (anyDuplicated(base$slug)) {
  dup <- base[duplicated(slug) | duplicated(slug, fromLast = TRUE)][order(slug)]
  stop(
    "Canonical species list produces duplicate slugs (needs disambiguation).\n",
    "First duplicates:\n",
    paste0(" - ", dup$species[1:min(10, nrow(dup))], collapse = "\n")
  )
}

# ==============================================================================
# Sensitive list (optional): flag sensitive taxa
# ==============================================================================

sens_abs <- find_by_stem(file.path(repo_root, sensitive_list_stem))
sens <- NULL
sens_slugs <- character()

if (!is.na(sens_abs) && file.exists(sens_abs)) {
  sens <- safe_read_csv_dt(sens_abs)
  
  if (!is.null(sens) && nrow(sens) > 0) {
    c_species <- pick_species_col(sens)
    
    spp_raw <- clean_text(sens[[c_species]])
    spp_raw <- spp_raw[!is.na(spp_raw) & nzchar(spp_raw)]
    
    spp_binom <- vapply(spp_raw, extract_binomial, FUN.VALUE = character(1))
    spp_binom <- spp_binom[!is.na(spp_binom) & nzchar(spp_binom)]
    
    sens_slugs <- unique(slugify_species(spp_binom))
  }
}

base[, sensitive_listed := slug %in% sens_slugs]

cat("[master_status] Canonical species rows: ", nrow(base), "\n", sep = "")
cat("[master_status] Sensitive species matched (canonical): ", sum(base$sensitive_listed, na.rm = TRUE), "\n", sep = "")
if (is.na(sens_abs) || !file.exists(sens_abs)) {
  cat("[master_status] NOTE: Sensitive list not found; sensitive_listed will be blank.\n")
  cat("[master_status] Expected stem: ", file.path(repo_root, sensitive_list_stem), "(.csv)\n", sep = "")
} else {
  cat("[master_status] Sensitive list used: ", sens_abs, "\n", sep = "")
}

# ==============================================================================
# Pipeline counts from runlogs (canonical-preferred)
# ==============================================================================

# ---- Stage 02 merge runlog (GBIF/NBN counts + merged totals) ------------------

merge_runlog <- file.path(processed_root, "02_merged", "_runlog_02_merged.csv")
if (!file.exists(merge_runlog)) merge_runlog <- file.path(processed_root, "01_merged", "_runlog_01_merged.csv")
if (!file.exists(merge_runlog)) {
  merge_runlog <- find_latest_by_basename_regex(
    processed_root,
    c("^_runlog_02_merged\\.csv$", "^_runlog_01_merged\\.csv$")
  )
}
merge_dt <- safe_read_csv_dt(merge_runlog)

merge_counts <- data.table(slug = base$slug)
merge_files <- data.table(
  slug = base$slug,
  gbif_file = NA_character_,
  nbn_file = NA_character_,
  gbif_exists = NA,
  nbn_exists = NA
)

if (!is.null(merge_dt) && nrow(merge_dt) > 0) {
  # Keep last row per slug (runlogs may contain multiple runs)
  if ("slug" %in% names(merge_dt)) {
    setorder(merge_dt, slug)
    merge_dt <- merge_dt[!duplicated(slug, fromLast = TRUE)]
  }
  
  pick <- function(dt, cands) {
    nm <- cands[cands %in% names(dt)]
    if (length(nm) == 0) return(NULL)
    nm[1]
  }
  
  c_slug <- pick(merge_dt, c("slug", "species_slug"))
  if (!is.null(c_slug) && c_slug != "slug") setnames(merge_dt, c_slug, "slug")
  
  c_gbif <- pick(merge_dt, c("n_gbif", "gbif_n", "n_gbif_after_strict_dedup"))
  c_nbn  <- pick(merge_dt, c("n_nbn", "nbn_n", "n_nbn_after_strict_dedup"))
  c_pre  <- pick(merge_dt, c("n_premerge", "n_total_after_strict_dedup", "n_total"))
  c_fin  <- pick(merge_dt, c("n_final", "n_after_merge", "n_merged"))
  
  merge_counts <- merge_dt[, .(
    slug = slug,
    gbif_records_clean_n = if (!is.null(c_gbif)) as_int(get(c_gbif)) else NA_integer_,
    nbn_records_clean_n  = if (!is.null(c_nbn))  as_int(get(c_nbn))  else NA_integer_,
    downloaded_total_before_cross_source_dedup_n = if (!is.null(c_pre)) as_int(get(c_pre)) else NA_integer_,
    merged_total_after_cross_source_dedup_n = if (!is.null(c_fin)) as_int(get(c_fin)) else NA_integer_
  )]
  
  # Fill pre-merge total if missing but sources exist
  merge_counts[
    is.na(downloaded_total_before_cross_source_dedup_n) &
      !is.na(gbif_records_clean_n) & !is.na(nbn_records_clean_n),
    downloaded_total_before_cross_source_dedup_n := gbif_records_clean_n + nbn_records_clean_n
  ]
  
  # Keep file paths + existence flags so we can backfill counts from disk if runlog counts are NA
  merge_files <- merge_dt[, .(
    slug = slug,
    gbif_file = if ("gbif_file" %in% names(merge_dt)) as.character(gbif_file) else NA_character_,
    nbn_file  = if ("nbn_file"  %in% names(merge_dt)) as.character(nbn_file)  else NA_character_,
    gbif_exists = if ("gbif_exists" %in% names(merge_dt)) as.logical(gbif_exists) else NA,
    nbn_exists  = if ("nbn_exists"  %in% names(merge_dt)) as.logical(nbn_exists)  else NA
  )]
  
  cat("[master_status] Merge runlog used: ", merge_runlog, "\n", sep = "")
  if ("slug" %in% names(merge_dt)) cat("[master_status] Merge runlog unique slugs: ", uniqueN(merge_dt$slug), "\n", sep = "")
} else {
  cat("[master_status] NOTE: Merge runlog not found/readable; merge counts will be NA.\n")
  cat("[master_status] Merge runlog resolved to: ", merge_runlog, "\n", sep = "")
}

# ---- Stage 04 filter runlog (policy filter in/out) ----------------------------

filter_runlog <- file.path(processed_root, "04_filtered", "_runlog_04_filtered.csv")
if (!file.exists(filter_runlog)) filter_runlog <- file.path(processed_root, "03_filtered", "_runlog_03_filtered.csv")
if (!file.exists(filter_runlog)) {
  filter_runlog <- find_latest_by_basename_regex(
    processed_root,
    c("^_runlog_04_filtered\\.csv$", "^_runlog_03_filtered\\.csv$")
  )
}
filter_dt <- safe_read_csv_dt(filter_runlog)

filter_counts <- data.table(slug = base$slug)

if (!is.null(filter_dt) && nrow(filter_dt) > 0) {
  if ("slug" %in% names(filter_dt)) {
    setorder(filter_dt, slug)
    filter_dt <- filter_dt[!duplicated(slug, fromLast = TRUE)]
  }
  
  pick <- function(dt, cands) {
    nm <- cands[cands %in% names(dt)]
    if (length(nm) == 0) return(NULL)
    nm[1]
  }
  
  c_slug <- pick(filter_dt, c("slug", "species_slug"))
  if (!is.null(c_slug) && c_slug != "slug") setnames(filter_dt, c_slug, "slug")
  
  c_in  <- pick(filter_dt, c("n_in", "n_before", "n_read", "n_before_policy_filter"))
  c_out <- pick(filter_dt, c("n_out", "n_after", "n_filtered", "n_after_policy_filter"))
  
  filter_counts <- filter_dt[, .(
    slug = slug,
    policy_filter_input_n = if (!is.null(c_in))  as_int(get(c_in))  else NA_integer_,
    filtered_total_after_policy_n = if (!is.null(c_out)) as_int(get(c_out)) else NA_integer_
  )]
  
  cat("[master_status] Filter runlog used: ", filter_runlog, "\n", sep = "")
  if ("slug" %in% names(filter_dt)) cat("[master_status] Filter runlog unique slugs: ", uniqueN(filter_dt$slug), "\n", sep = "")
} else {
  cat("[master_status] NOTE: Filter runlog not found/readable; filter counts will be NA.\n")
  cat("[master_status] Filter runlog resolved to: ", filter_runlog, "\n", sep = "")
}

# ---- Grid summary (occupied cells + GEE points) -------------------------------

grid_summary <- find_latest_by_basename_regex(processed_root, c("^_summary_grid\\.csv$"))
grid_dt <- safe_read_csv_dt(grid_summary)

grid_counts <- data.table(slug = base$slug)

if (!is.null(grid_dt) && nrow(grid_dt) > 0) {
  pick <- function(dt, cands) {
    nm <- cands[cands %in% names(dt)]
    if (length(nm) == 0) return(NULL)
    nm[1]
  }
  
  c_slug <- pick(grid_dt, c("slug", "species_slug"))
  if (!is.null(c_slug) && c_slug != "slug") setnames(grid_dt, c_slug, "slug")
  
  c_cells <- pick(grid_dt, c("n_presence_cells", "presence_cells", "n_gridded_presence_cells"))
  c_pts   <- pick(grid_dt, c("n_points_for_gee", "n_gridded_points_for_gee", "gee_points", "n_gee_points"))
  
  grid_counts <- grid_dt[, .(
    slug = slug,
    gridded_occupied_cells_n = if (!is.null(c_cells)) as_int(get(c_cells)) else NA_integer_,
    gee_points_n = if (!is.null(c_pts)) as_int(get(c_pts)) else NA_integer_
  )]
  
  cat("[master_status] Grid summary used: ", grid_summary, "\n", sep = "")
} else {
  cat("[master_status] NOTE: Grid summary not found/readable; grid counts will be NA.\n")
}

# ==============================================================================
# Uncertainty / obscuring audit (pre-filter; Stage 02 merged parquet)
# ==============================================================================

merged_path_for <- function(slug) {
  file.path(processed_root, "02_merged", slug, paste0("occ_", slug, "__merged.parquet"))
}

audit_uncertainty_one <- function(slug, is_sensitive) {
  p <- merged_path_for(slug)
  
  if (!file.exists(p)) {
    return(data.table(
      slug = slug,
      merged_records_pre_filter_n = NA_integer_,
      uncertainty_threshold_used_m = as.integer(default_max_uncertainty_m),
      uncertainty_non_missing_pct_pre_filter = NA_real_,
      uncertainty_mode_m_pre_filter = NA_real_,
      would_drop_if_uncertainty_gt_threshold_n = NA_integer_,
      would_drop_if_uncertainty_gt_threshold_pct_of_non_missing = NA_real_,
      uk_records_pre_filter_n = NA_integer_,
      uk_uncertainty_mode_m_pre_filter = NA_real_
    ))
  }
  
  cols <- c("coordinateUncertaintyInMeters")
  if (isTRUE(compute_uk_stats_for_sensitive) && isTRUE(is_sensitive)) cols <- c(cols, "country")
  
  tab <- read_parquet_cols(p, cols)
  if (is.null(tab) || !"coordinateUncertaintyInMeters" %in% names(tab)) {
    return(data.table(
      slug = slug,
      merged_records_pre_filter_n = NA_integer_,
      uncertainty_threshold_used_m = as.integer(default_max_uncertainty_m),
      uncertainty_non_missing_pct_pre_filter = NA_real_,
      uncertainty_mode_m_pre_filter = NA_real_,
      would_drop_if_uncertainty_gt_threshold_n = NA_integer_,
      would_drop_if_uncertainty_gt_threshold_pct_of_non_missing = NA_real_,
      uk_records_pre_filter_n = NA_integer_,
      uk_uncertainty_mode_m_pre_filter = NA_real_
    ))
  }
  
  n_total <- nrow(tab)
  
  u <- suppressWarnings(as.numeric(tab$coordinateUncertaintyInMeters))
  is_ok <- is.finite(u)
  n_ok <- sum(is_ok)
  pct_ok <- if (n_total > 0) round(100 * n_ok / n_total, 1) else NA_real_
  
  u_ok <- u[is_ok]
  u_mode <- if (length(u_ok) > 0) mode_numeric(u_ok) else NA_real_
  
  # uncertainty-only keep/drop (mirrors the uncertainty rule)
  keep_unc <- rep(TRUE, n_total)
  keep_unc[is_ok] <- u[is_ok] <= default_max_uncertainty_m
  if (identical(uncertainty_missing_action, "drop")) keep_unc[!is_ok] <- FALSE
  
  n_drop_unc <- as.integer(n_total - sum(keep_unc))
  
  # % dropped among records with known (numeric) uncertainty
  n_drop_known <- if (length(u_ok) > 0) as.integer(sum(u_ok > default_max_uncertainty_m)) else 0L
  pct_drop_known <- if (length(u_ok) > 0) round(100 * n_drop_known / length(u_ok), 1) else NA_real_
  
  # UK-only mode (sensitive only)
  uk_n <- NA_integer_
  uk_mode <- NA_real_
  
  if (isTRUE(compute_uk_stats_for_sensitive) && isTRUE(is_sensitive) && ("country" %in% names(tab))) {
    is_uk <- tab$country %in% uk_country_values
    uk_n <- as.integer(sum(is_uk, na.rm = TRUE))
    
    if (uk_n > 0) {
      u_uk <- suppressWarnings(as.numeric(tab$coordinateUncertaintyInMeters[is_uk]))
      u_uk_ok <- u_uk[is.finite(u_uk)]
      if (length(u_uk_ok) > 0) uk_mode <- mode_numeric(u_uk_ok)
    }
  }
  
  data.table(
    slug = slug,
    merged_records_pre_filter_n = as.integer(n_total),
    uncertainty_threshold_used_m = as.integer(default_max_uncertainty_m),
    uncertainty_non_missing_pct_pre_filter = pct_ok,
    uncertainty_mode_m_pre_filter = u_mode,
    would_drop_if_uncertainty_gt_threshold_n = n_drop_unc,
    would_drop_if_uncertainty_gt_threshold_pct_of_non_missing = pct_drop_known,
    uk_records_pre_filter_n = uk_n,
    uk_uncertainty_mode_m_pre_filter = uk_mode
  )
}

cat("[master_status] Auditing uncertainty from Stage 02 merged parquet...\n")
t0 <- Sys.time()

u_audit <- rbindlist(lapply(seq_len(nrow(base)), function(i) {
  if (i %% 50 == 0) cat("[master_status] ", i, "/", nrow(base), "\n", sep = "")
  audit_uncertainty_one(base$slug[i], base$sensitive_listed[i])
}), fill = TRUE)

cat("[master_status] Uncertainty audit runtime (sec): ",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n", sep = "")

# ==============================================================================
# Assemble master table
# ==============================================================================

st <- copy(base)

st <- merge(st, merge_counts, by = "slug", all.x = TRUE)
st <- merge(st, merge_files,  by = "slug", all.x = TRUE)
st <- merge(st, filter_counts, by = "slug", all.x = TRUE)
st <- merge(st, grid_counts, by = "slug", all.x = TRUE)
st <- merge(st, u_audit, by = "slug", all.x = TRUE)

# ------------------------------------------------------------------------------
# Backfill merge-stage counts from disk (robust to partial reruns / NA runlogs)
# ------------------------------------------------------------------------------

if (isTRUE(backfill_counts_from_files)) {
  
  resolve_to_abs <- function(p) {
    p <- as.character(p)
    if (is.na(p) || !nzchar(p)) return(NA_character_)
    if (file.exists(p)) return(p)
    p2 <- file.path(repo_root, p)
    if (file.exists(p2)) return(p2)
    NA_character_
  }
  
  count_lines_fast <- function(path) {
    if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_integer_)
    sz <- file.info(path)$size
    if (is.na(sz) || sz == 0) return(0L)
    
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    
    n_nl <- 0L
    repeat {
      buf <- readBin(con, what = "raw", n = 1024L * 1024L)
      if (length(buf) == 0) break
      n_nl <- n_nl + sum(buf == as.raw(10))
    }
    
    # If last line has no trailing newline, add one
    last_is_nl <- FALSE
    con2 <- file(path, open = "rb")
    on.exit(close(con2), add = TRUE)
    seek(con2, where = max(0, sz - 1), origin = "start")
    last <- readBin(con2, what = "raw", n = 1)
    last_is_nl <- length(last) == 1 && identical(last, as.raw(10))
    
    as.integer(n_nl + ifelse(last_is_nl, 0L, 1L))
  }
  
  count_rows_text_with_header <- function(path) {
    n <- count_lines_fast(path)
    if (is.na(n)) return(NA_integer_)
    as.integer(max(0L, n - 1L))
  }
  
  count_rows_any <- function(path) {
    if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_integer_)
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("parquet")) {
      return(as.integer(tryCatch(arrow::read_parquet(path, col_select = 1)$length(), error = function(e) NA_integer_)))
    }
    if (ext %in% c("csv", "tsv", "txt")) return(count_rows_text_with_header(path))
    count_rows_text_with_header(path)
  }
  
  # Resolve paths once
  st[, gbif_file_abs := vapply(gbif_file, resolve_to_abs, FUN.VALUE = character(1))]
  st[, nbn_file_abs  := vapply(nbn_file,  resolve_to_abs, FUN.VALUE = character(1))]
  
  # Fill GBIF counts
  need_gbif <- which(is.na(st$gbif_records_clean_n) & !is.na(st$gbif_file_abs))
  if (length(need_gbif) > 0) {
    st[need_gbif, gbif_records_clean_n := vapply(gbif_file_abs, count_rows_any, FUN.VALUE = integer(1))]
  }
  st[is.na(gbif_records_clean_n) & isFALSE(gbif_exists), gbif_records_clean_n := 0L]
  
  # Fill NBN counts
  need_nbn <- which(is.na(st$nbn_records_clean_n) & !is.na(st$nbn_file_abs))
  if (length(need_nbn) > 0) {
    st[need_nbn, nbn_records_clean_n := vapply(nbn_file_abs, count_rows_any, FUN.VALUE = integer(1))]
  }
  st[is.na(nbn_records_clean_n) & isFALSE(nbn_exists), nbn_records_clean_n := 0L]
  
  # Fill downloaded total if missing
  st[
    is.na(downloaded_total_before_cross_source_dedup_n) &
      !is.na(gbif_records_clean_n) & !is.na(nbn_records_clean_n),
    downloaded_total_before_cross_source_dedup_n := gbif_records_clean_n + nbn_records_clean_n
  ]
  
  # Fill merged total if missing using the Stage 02 merged parquet row count (already computed in u_audit)
  st[
    is.na(merged_total_after_cross_source_dedup_n) & !is.na(merged_records_pre_filter_n),
    merged_total_after_cross_source_dedup_n := merged_records_pre_filter_n
  ]
  
  # Clean up internal helper columns
  st[, c("gbif_file_abs", "nbn_file_abs") := NULL]
  
  cat("[master_status] Backfill complete: non-missing gbif_records_clean_n = ",
      sum(!is.na(st$gbif_records_clean_n)), "/", nrow(st), "\n", sep = "")
  cat("[master_status] Backfill complete: non-missing merged_total_after_cross_source_dedup_n = ",
      sum(!is.na(st$merged_total_after_cross_source_dedup_n)), "/", nrow(st), "\n", sep = "")
}

# Derived metrics: cross-source de-duplication
st[, duplicates_removed_cross_source_n := ifelse(
  !is.na(downloaded_total_before_cross_source_dedup_n) & !is.na(merged_total_after_cross_source_dedup_n),
  downloaded_total_before_cross_source_dedup_n - merged_total_after_cross_source_dedup_n,
  NA_integer_
)]
st[, duplicates_removed_cross_source_pct := pct1(
  duplicates_removed_cross_source_n,
  downloaded_total_before_cross_source_dedup_n
)]

# Derived metrics: policy filtering (overall)
st[, dropped_by_policy_total_n := ifelse(
  !is.na(policy_filter_input_n) & !is.na(filtered_total_after_policy_n),
  policy_filter_input_n - filtered_total_after_policy_n,
  NA_integer_
)]
st[, dropped_by_policy_total_pct := pct1(dropped_by_policy_total_n, policy_filter_input_n)]

# Extra drops beyond uncertainty-only, using best available merged total:
# prefer merged_total_after_cross_source_dedup_n; otherwise fall back to merged_records_pre_filter_n
st[, merged_total_for_uncertainty_math_n := merged_total_after_cross_source_dedup_n]
st[is.na(merged_total_for_uncertainty_math_n) & !is.na(merged_records_pre_filter_n),
   merged_total_for_uncertainty_math_n := merged_records_pre_filter_n]

st[, uncertainty_kept_only_n := ifelse(
  !is.na(merged_total_for_uncertainty_math_n) & !is.na(would_drop_if_uncertainty_gt_threshold_n),
  merged_total_for_uncertainty_math_n - would_drop_if_uncertainty_gt_threshold_n,
  NA_integer_
)]

st[, extra_drops_from_other_policy_rules_n := ifelse(
  !is.na(uncertainty_kept_only_n) & !is.na(filtered_total_after_policy_n),
  uncertainty_kept_only_n - filtered_total_after_policy_n,
  NA_integer_
)]
st[, extra_drops_from_other_policy_rules_pct_of_uncertainty_kept := pct1(
  extra_drops_from_other_policy_rules_n,
  uncertainty_kept_only_n
)]

# Stage reached (no diagnosis/status fields)
st[, stage_reached := fifelse(
  !is.na(gridded_occupied_cells_n), "Gridded",
  fifelse(!is.na(filtered_total_after_policy_n), "Filtered",
          fifelse(!is.na(merged_total_after_cross_source_dedup_n) | !is.na(merged_records_pre_filter_n), "Merged",
                  fifelse(!is.na(downloaded_total_before_cross_source_dedup_n), "Downloaded", "Unknown")
          )
  )
)]

# Keep simple GEE placeholders for manual tracking
st[, gee_assets_exist := ""]
st[, gee_model_run := ""]

# ==============================================================================
# Final output table (sharing-friendly)
# ==============================================================================

out <- st[, .(
  species,
  slug,
  sensitive_listed = ifelse(slug %in% sens_slugs, "Yes", "No"),
  
  stage_reached,
  
  # Data volumes
  gbif_records_clean_n,
  nbn_records_clean_n,
  downloaded_total_before_cross_source_dedup_n,
  
  merged_total_after_cross_source_dedup_n,
  duplicates_removed_cross_source_n,
  duplicates_removed_cross_source_pct,
  
  filtered_total_after_policy_n,
  dropped_by_policy_total_n,
  dropped_by_policy_total_pct,
  
  gridded_occupied_cells_n,
  gee_points_n,
  
  # Uncertainty / obscuring impact (pre-filter, uncertainty-only)
  uncertainty_threshold_used_m,
  uncertainty_non_missing_pct_pre_filter,
  uncertainty_mode_m_pre_filter,
  would_drop_if_uncertainty_gt_threshold_n,
  would_drop_if_uncertainty_gt_threshold_pct_of_non_missing,
  extra_drops_from_other_policy_rules_n,
  extra_drops_from_other_policy_rules_pct_of_uncertainty_kept,
  
  # UK-only (sensitive species only; blank otherwise)
  uk_records_pre_filter_n,
  uk_uncertainty_mode_m_pre_filter,
  
  # Placeholders for later tracking
  gee_assets_exist,
  gee_model_run
)]

# Sort: sensitive first, then biggest uncertainty-rule drop, then lowest gridded cells
out[, sens_rank := ifelse(sensitive_listed == "Yes", 1L, 0L)]
setorder(out, -sens_rank, -would_drop_if_uncertainty_gt_threshold_n, gridded_occupied_cells_n, na.last = TRUE)
out[, sens_rank := NULL]

# ==============================================================================
# Write CSV + README
# ==============================================================================

out_abs <- file.path(repo_root, out_csv)
dir.create(dirname(out_abs), recursive = TRUE, showWarnings = FALSE)
fwrite(out, out_abs, na = "")

if (isTRUE(write_readme)) {
  readme_path <- sub("\\.csv$", "__README.txt", out_abs)
  
  # Per-column glossary: must cover every header in the CSV
  defs <- list(
    species = "Species scientific name (binomial) from the canonical species list in data/_meta/species_list_binomial.csv (no header; one binomial per line).",
    slug = "Lowercase, underscore-separated species identifier derived from the binomial; used for folder/file naming.",
    sensitive_listed = "Yes/No flag: species appears in the project’s sensitive-species list (if found).",
    
    stage_reached = "Furthest completed stage inferred from whether key outputs exist for the species: Downloaded → Merged → Filtered → Gridded.",
    
    gbif_records_clean_n = "Number of GBIF records after GBIF-specific cleaning and strict de-duplication (pre-merge). Read from merge runlog when available; otherwise backfilled by counting rows in the GBIF clean file if available.",
    nbn_records_clean_n = "Number of NBN records after NBN-specific cleaning and strict de-duplication (pre-merge). Read from merge runlog when available; otherwise backfilled by counting rows in the NBN clean file if available.",
    downloaded_total_before_cross_source_dedup_n = "gbif_records_clean_n + nbn_records_clean_n (pre-merge total), if available.",
    
    merged_total_after_cross_source_dedup_n = "Number of records after combining GBIF and NBN and removing duplicates between sources. Read from merge runlog when available; otherwise backfilled using the Stage 02 merged parquet row count (pre-filter).",
    duplicates_removed_cross_source_n = "Records removed during cross-source de-duplication: downloaded_total_before_cross_source_dedup_n - merged_total_after_cross_source_dedup_n.",
    duplicates_removed_cross_source_pct = "Percentage removed during cross-source de-duplication: duplicates_removed_cross_source_n / downloaded_total_before_cross_source_dedup_n * 100.",
    
    filtered_total_after_policy_n = "Number of records remaining after the Stage 04 policy filter (all rules), from the filter runlog (if available).",
    dropped_by_policy_total_n = "Records removed by the Stage 04 policy filter overall: policy_filter_input_n - filtered_total_after_policy_n (policy_filter_input_n comes from the filter runlog if available).",
    dropped_by_policy_total_pct = "Percentage removed by the Stage 04 policy filter overall: dropped_by_policy_total_n / policy_filter_input_n * 100.",
    
    gridded_occupied_cells_n = "Number of occupied grid cells after gridding (unique presence cells), from the latest _summary_grid.csv found under data/processed (if available).",
    gee_points_n = "Number of point features exported for modelling, from the latest _summary_grid.csv (if available). Typically one point per occupied cell using the first-location-in-cell rule; can be lower if points are dropped by masking or missing predictor coverage.",
    
    uncertainty_threshold_used_m = "Uncertainty threshold (meters) used for the uncertainty-only audit; intended to match the Stage 04 policy setting (e.g., 1000 m).",
    uncertainty_non_missing_pct_pre_filter = "Percentage of pre-filter merged records with a numeric coordinateUncertaintyInMeters value.",
    uncertainty_mode_m_pre_filter = "Most common (modal) numeric coordinateUncertaintyInMeters value in pre-filter merged records. A spike near ~7071 m can indicate systematic 10 km-square generalisation (5 km * sqrt(2)).",
    would_drop_if_uncertainty_gt_threshold_n = "Number of pre-filter merged records that would be dropped by the uncertainty rule alone (applying only the threshold; ignoring other policy rules).",
    would_drop_if_uncertainty_gt_threshold_pct_of_non_missing = "Of records with a numeric uncertainty value, the percentage exceeding the threshold.",
    extra_drops_from_other_policy_rules_n = "Additional records dropped by other policy rules after accounting for uncertainty-only dropping: (uncertainty-kept-only) - filtered_total_after_policy_n.",
    extra_drops_from_other_policy_rules_pct_of_uncertainty_kept = "Percentage of uncertainty-kept-only records that are further dropped by other policy rules.",
    
    uk_records_pre_filter_n = "Pre-filter merged record count where country is one of the UK country values (populated only for sensitive-listed species).",
    uk_uncertainty_mode_m_pre_filter = "Modal uncertainty value for UK-only pre-filter merged records (populated only for sensitive-listed species).",
    
    gee_assets_exist = "Placeholder (manual): whether expected GEE assets exist for this species.",
    gee_model_run = "Placeholder (manual): whether a model has been run for this species."
  )
  
  gloss <- character()
  gloss <- c(gloss, "InfluentialSpecies — Master species pipeline status (with uncertainty/obscuring audit)")
  gloss <- c(gloss, "===========================================================================")
  gloss <- c(gloss, "")
  gloss <- c(gloss, "Row set")
  gloss <- c(gloss, "-------")
  gloss <- c(gloss, "Rows are defined ONLY by: data/_meta/species_list_binomial.csv (no header; one binomial per line).")
  gloss <- c(gloss, "")
  gloss <- c(gloss, "Key field")
  gloss <- c(gloss, "---------")
  gloss <- c(gloss, "coordinateUncertaintyInMeters is treated as an uncertainty radius (meters) around the reported coordinates.")
  gloss <- c(gloss, "")
  gloss <- c(gloss, "10 km-square generalisation signal (example)")
  gloss <- c(gloss, "------------------------------------------")
  gloss <- c(gloss, "A common uncertainty value near 7071 m can reflect 10 km-square generalisation.")
  gloss <- c(gloss, "If coordinates are generalised to a 10 km square, the half-side length is 5 km.")
  gloss <- c(gloss, "The distance from the square centre to a corner is: sqrt(5^2 + 5^2) km = 5 * sqrt(2) km ≈ 7.071 km = 7071 m.")
  gloss <- c(gloss, "")
  gloss <- c(gloss, "Uncertainty-only audit behaviour")
  gloss <- c(gloss, "-------------------------------")
  gloss <- c(gloss, paste0("uncertainty_threshold_used_m = ", default_max_uncertainty_m, " (meters)"))
  gloss <- c(gloss, paste0("uncertainty_missing_action  = '", uncertainty_missing_action, "'"))
  gloss <- c(gloss, "")
  gloss <- c(gloss, "Column glossary")
  gloss <- c(gloss, "--------------")
  gloss <- c(gloss, "")
  
  missing_defs <- character()
  for (nm in names(out)) {
    if (!nm %in% names(defs)) {
      missing_defs <- c(missing_defs, nm)
      gloss <- c(gloss, nm, "  (No definition available in this README for this column.)", "")
    } else {
      gloss <- c(gloss, nm, paste0("  ", defs[[nm]]), "")
    }
  }
  
  if (length(missing_defs) > 0) {
    cat("[master_status] WARNING: README missing definitions for columns:\n")
    cat(paste0(" - ", missing_defs, collapse = "\n"), "\n", sep = "")
  }
  
  gloss <- c(gloss, "Files written by this script")
  gloss <- c(gloss, "----------------------------")
  gloss <- c(gloss, paste0("  - ", out_csv))
  gloss <- c(gloss, paste0("  - ", sub("\\.csv$", "__README.txt", out_csv)))
  gloss <- c(gloss, "")
  
  writeLines(gloss, readme_path)
}

cat("\n============================================================\n")
cat("Master pipeline status written\n")
cat("============================================================\n")
cat("Output: ", out_abs, "\n", sep = "")
if (isTRUE(write_readme)) cat("README: ", sub("\\.csv$", "__README.txt", out_abs), "\n", sep = "")
cat("Species total: ", nrow(out), "\n", sep = "")
cat("Sensitive matched: ", sum(out$sensitive_listed == "Yes", na.rm = TRUE), "\n", sep = "")
cat("============================================================\n\n")