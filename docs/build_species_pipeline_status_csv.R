# docs/build_species_pipeline_status_csv.R --------------------------------------
#
# Purpose
#   Build a per-species status table (CSV) covering the
#   pipeline up to (and including) gridding for GEE.
#
# What it tries to capture (one row per species)
#   - Species list (binomial) + slug
#   - Stage 02: merge runlog (GBIF/NBN presence + counts; after cross-source de-dup)
#   - Stage 03: QC flagging status + input count
#   - Stage 04: policy filtering status + n_in / n_out + policy_id
#   - Stage 05: gridding summary (presence cells + GEE points)
#   - Placeholders for later: GEE assets/model status, sensitive-species losses, etc.
#
# Inputs
#   data/_meta/species_list_binomial.csv
#   data/processed/* runlogs (merge / qc / filter), plus:
#   data/processed/**/_summary_grid.csv
#
# Output
#   docs/derived/species_pipeline_status.csv
#
# Restart safety
#   - Stage 02 merge runlog may contain "skipped_cached" rows with NA counts.
#   - This script back-fills missing counts by counting rows in the referenced files.
#   - To make that resumable, per-file counts are cached to:
#       docs/derived/species_pipeline_status__count_cache.csv
#     If you stop and re-run, already-counted files (unchanged size/mtime) are skipped.
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Repo root (robust when sourced from docs/ or scripts/) -------------------
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

find_repo_root <- function(start_dir) {
  markers <- c(".git", "R", "data", "InfluentialSpecies.Rproj", "DESCRIPTION")
  d <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:15) {
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

# ---- Output folder + cache ----------------------------------------------------
out_dir <- file.path(repo_root, "docs", "derived")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_path <- file.path(out_dir, "species_pipeline_status.csv")

# Append-only cache of per-file row counts (CSV). Safe to delete if needed.
cache_path <- file.path(out_dir, "species_pipeline_status__count_cache.csv")

# ---- Small helpers ------------------------------------------------------------
slugify_species <- function(x) {
  s <- gsub("[^a-z0-9]+", "_", tolower(as.character(x)))
  gsub("^_+|_+$", "", s)
}

yesno <- function(x) {
  ifelse(is.na(x), "", ifelse(isTRUE(x), "Yes", "No"))
}

first_existing <- function(paths) {
  paths <- unique(paths[nzchar(paths)])
  for (p in paths) if (file.exists(p)) return(p)
  NA_character_
}

pick_col <- function(dt, candidates) {
  nm <- candidates[candidates %in% names(dt)]
  if (length(nm) == 0) return(NULL)
  nm[[1]]
}

safe_read_csv_dt <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  dt <- tryCatch(
    fread(path, sep = ",", header = TRUE, fill = TRUE, quote = "\"", na.strings = c("", "NA")),
    error = function(e) NULL
  )
  if (is.null(dt)) return(NULL)
  
  hdr <- tryCatch(strsplit(readLines(path, n = 1, warn = FALSE), ",", fixed = TRUE)[[1]], error = function(e) NULL)
  if (!is.null(hdr) && ncol(dt) > length(hdr)) {
    setnames(dt, c(hdr, paste0("extra_col_", seq_len(ncol(dt) - length(hdr)))))
  }
  dt
}

latest_by_slug <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  
  slug_col <- pick_col(dt, c("slug", "species_slug"))
  if (is.null(slug_col)) return(dt[0])
  
  ts_col <- pick_col(dt, c("timestamp_utc", "timestamp", "time_utc"))
  dt[, .row_id__ := .I]
  
  if (!is.null(ts_col)) {
    dt[, .ts__ := suppressWarnings(as.POSIXct(get(ts_col), format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))]
    if (all(is.na(dt$.ts__))) dt[, .ts__ := suppressWarnings(as.POSIXct(get(ts_col), tz = "UTC"))]
    setorderv(dt, c(".ts__", ".row_id__"), c(1, 1), na.last = TRUE)
  } else {
    setorder(dt, .row_id__)
  }
  
  out <- dt[, .SD[.N], by = slug_col]
  setnames(out, slug_col, "slug")
  
  drop_cols <- intersect(c(".row_id__", ".ts__"), names(out))
  if (length(drop_cols) > 0) out[, (drop_cols) := NULL]
  
  out
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

fix_path <- function(p) {
  p <- trimws(as.character(p))
  p[p == ""] <- NA_character_
  p
}

# ---- Count cache helpers (restart-safe) --------------------------------------
load_count_cache <- function(path) {
  if (!file.exists(path)) {
    return(data.table(path = character(), size = numeric(), mtime = character(), n_rows = integer()))
  }
  cc <- tryCatch(fread(path, fill = TRUE), error = function(e) NULL)
  if (is.null(cc) || nrow(cc) == 0) {
    return(data.table(path = character(), size = numeric(), mtime = character(), n_rows = integer()))
  }
  if (!all(c("path", "size", "mtime", "n_rows") %in% names(cc))) {
    return(data.table(path = character(), size = numeric(), mtime = character(), n_rows = integer()))
  }
  cc[, path := as.character(path)]
  cc[, size := suppressWarnings(as.numeric(size))]
  cc[, mtime := as.character(mtime)]
  cc[, n_rows := suppressWarnings(as.integer(n_rows))]
  cc[!is.na(path) & nzchar(path)]
}

append_count_cache <- function(path, rec) {
  # rec: data.table with columns path,size,mtime,n_rows
  # Ensure header appears once, then append data only.
  write_header <- !file.exists(path)
  fwrite(rec, path, append = !write_header, col.names = write_header)
}

count_csv_rows_fast <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_integer_)
  
  con <- tryCatch(file(path, open = "rb"), error = function(e) NULL)
  if (is.null(con)) return(NA_integer_)
  
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  
  buf_size <- 8L * 1024L * 1024L
  n_nl <- 0L
  
  repeat {
    raw <- tryCatch(readBin(con, what = "raw", n = buf_size), error = function(e) NULL)
    if (is.null(raw) || length(raw) == 0) break
    n_nl <- n_nl + sum(raw == as.raw(10))  # '\n'
  }
  
  as.integer(max(0L, n_nl - 1L))  # subtract header line
}

count_rows_table_file <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_integer_)
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) return(NA_integer_)
    tab <- tryCatch(arrow::read_parquet(path), error = function(e) NULL)
    if (is.null(tab)) return(NA_integer_)
    return(as.integer(nrow(tab)))
  }
  
  if (ext == "rds") {
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(obj)) return(NA_integer_)
    return(as.integer(NROW(obj)))
  }
  
  NA_integer_
}

cached_count <- function(cc, path_in) {
  # Always return a scalar integer NA or count.
  if (is.null(path_in) || length(path_in) != 1L) return(NA_integer_)
  
  if (is.na(path_in) || !nzchar(path_in) || !file.exists(path_in)) return(NA_integer_)
  if (is.null(cc) || !inherits(cc, "data.table") || nrow(cc) == 0) return(NA_integer_)
  if (!all(c("path", "size", "mtime", "n_rows") %in% names(cc))) return(NA_integer_)
  
  info <- file.info(path_in)
  if (is.na(info$size) || is.na(info$mtime)) return(NA_integer_)
  
  hit <- cc[
    path == path_in &
      size == as.numeric(info$size) &
      mtime == as.character(info$mtime)
  ]
  
  if (nrow(hit) == 0) return(NA_integer_)
  
  v <- hit$n_rows[.N]
  if (is.null(v) || length(v) != 1L) return(NA_integer_)
  as.integer(v)
}

set_cached_count <- function(cc, file_path, n_rows, cache_path) {
  if (is.null(file_path) || length(file_path) != 1L) return(cc)
  if (is.na(file_path) || !nzchar(file_path) || !file.exists(file_path)) return(cc)
  
  info <- file.info(file_path)
  
  rec <- data.table(
    path = as.character(file_path),
    size = as.numeric(info$size),
    mtime = as.character(info$mtime),
    n_rows = as.integer(n_rows)
  )
  
  # In-memory update: drop older entries for same file, then append
  if (is.null(cc) || !inherits(cc, "data.table")) {
    cc <- data.table(path = character(), size = numeric(), mtime = character(), n_rows = integer())
  }
  if (!all(c("path", "size", "mtime", "n_rows") %in% names(cc))) {
    cc <- data.table(path = character(), size = numeric(), mtime = character(), n_rows = integer())
  }
  
  cc <- cc[path != rec$path]
  cc <- rbind(cc, rec, fill = TRUE)
  
  # On-disk append
  append_count_cache(cache_path, rec)
  
  cc
}

count_csv_rows_cached <- function(path, cc, cache_path) {
  n0 <- cached_count(cc, path)
  if (is.null(n0) || length(n0) != 1L) n0 <- NA_integer_
  if (!is.na(n0)) return(list(n = n0, cc = cc, used_cache = TRUE))
  
  n <- count_csv_rows_fast(path)
  cc2 <- set_cached_count(cc, path, n, cache_path)
  list(n = n, cc = cc2, used_cache = FALSE)
}

count_table_rows_cached <- function(path, cc, cache_path) {
  n0 <- cached_count(cc, path)
  if (is.null(n0) || length(n0) != 1L) n0 <- NA_integer_
  if (!is.na(n0)) return(list(n = n0, cc = cc, used_cache = TRUE))
  
  n <- count_rows_table_file(path)
  cc2 <- set_cached_count(cc, path, n, cache_path)
  list(n = n, cc = cc2, used_cache = FALSE)
}

# ---- Species list (canonical) -------------------------------------------------
meta_species_csv <- file.path(repo_root, "data", "_meta", "species_list_binomial.csv")
if (!file.exists(meta_species_csv)) stop("Can't find: ", meta_species_csv)

binom <- tryCatch(
  read.csv(meta_species_csv, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) NULL
)
if (is.null(binom) || nrow(binom) == 0) stop("Couldn't read species list: ", meta_species_csv)

species_names <- trimws(as.character(binom[[1]]))
species_names <- species_names[!is.na(species_names) & nzchar(species_names)]
species_names <- species_names[tolower(species_names) != "binomial"]
species_names <- species_names[!duplicated(species_names)]
if (length(species_names) < 50) stop("Species list looks unexpectedly short (", length(species_names), "). Check: ", meta_species_csv)

base <- data.table(
  species_binomial = species_names,
  slug = slugify_species(species_names)
)

# ---- Locate runlogs -----------------------------------------------------------
merge_runlog <- first_existing(c(
  file.path(processed_root, "02_merged", "_runlog_02_merged.csv"),
  file.path(processed_root, "01_merged", "_runlog_01_merged.csv")
))
if (is.na(merge_runlog)) merge_runlog <- find_latest_by_basename_regex(processed_root, c("^_runlog_02_merged\\.csv$", "^_runlog_01_merged\\.csv$"))

qc_runlog <- first_existing(c(
  file.path(processed_root, "03_qc_flagged", "_runlog_03_qc_flagged.csv"),
  file.path(processed_root, "02_qc_flagged", "_runlog_02_qc_flagged.csv")
))
if (is.na(qc_runlog)) qc_runlog <- find_latest_by_basename_regex(processed_root, c("^_runlog_03_qc_flagged\\.csv$", "^_runlog_02_qc_flagged\\.csv$"))

filter_runlog <- first_existing(c(
  file.path(processed_root, "04_filtered", "_runlog_04_filtered.csv"),
  file.path(processed_root, "03_filtered", "_runlog_03_filtered.csv")
))
if (is.na(filter_runlog)) filter_runlog <- find_latest_by_basename_regex(processed_root, c("^_runlog_04_filtered\\.csv$", "^_runlog_03_filtered\\.csv$"))

grid_runlog <- find_latest_by_basename_regex(processed_root, c("^_runlog_05_grid\\.csv$", "^_runlog_04_grid\\.csv$"))
grid_summary <- find_latest_by_basename_regex(processed_root, c("^_summary_grid\\.csv$"))

# ---- Load / init cache --------------------------------------------------------
cc <- load_count_cache(cache_path)
if (nrow(cc) > 0) {
  cat("[status_table] Loaded count cache rows: ", nrow(cc), "\n", sep = "")
} else {
  cat("[status_table] No count cache yet (will create): ", cache_path, "\n", sep = "")
}

# ---- Stage 02 merge -----------------------------------------------------------
merge_dt <- safe_read_csv_dt(merge_runlog)
merge_latest <- latest_by_slug(merge_dt)

merge_std <- NULL
if (!is.null(merge_latest) && nrow(merge_latest) > 0) {
  
  merge_latest[, gbif_file := fix_path(gbif_file)]
  merge_latest[, nbn_file  := fix_path(nbn_file)]
  merge_latest[, out_file  := fix_path(out_file)]
  
  merge_std <- merge_latest[, .(
    slug = slug,
    
    raw_group_dir_used = as.character(group_dir),
    
    gbif_downloaded = as.logical(gbif_exists),
    n_gbif_after_strict_dedup = suppressWarnings(as.integer(n_gbif)),
    gbif_file = as.character(gbif_file),
    
    nbn_downloaded = as.logical(nbn_exists),
    n_nbn_after_strict_dedup = suppressWarnings(as.integer(n_nbn)),
    nbn_file = as.character(nbn_file),
    
    n_total_after_strict_dedup = suppressWarnings(as.integer(n_premerge)),
    n_dropped_between_sources = suppressWarnings(as.integer(n_dropped_strict)),
    n_after_merge = suppressWarnings(as.integer(n_final)),
    
    out_file = as.character(out_file),
    stage02_merge_status = as.character(status),
    stage02_merge_note = as.character(note)
  )]
  
  # Back-fill missing counts for cached skips (resumable via cache).
  needs_gbif <- which(is.na(merge_std$n_gbif_after_strict_dedup) & !is.na(merge_std$gbif_file) & file.exists(merge_std$gbif_file))
  needs_nbn  <- which(is.na(merge_std$n_nbn_after_strict_dedup)  & !is.na(merge_std$nbn_file)  & file.exists(merge_std$nbn_file))
  needs_out  <- which(is.na(merge_std$n_after_merge)             & !is.na(merge_std$out_file)  & file.exists(merge_std$out_file))
  
  cat("\n[status_table] Stage 02 back-fill needed (files exist, counts missing):\n")
  cat(" - GBIF clean CSVs: ", length(needs_gbif), "\n", sep = "")
  cat(" - NBN  clean CSVs: ", length(needs_nbn),  "\n", sep = "")
  cat(" - Merged outputs:  ", length(needs_out),  "\n\n", sep = "")
  
  if (length(needs_gbif) > 0) {
    for (k in seq_along(needs_gbif)) {
      i <- needs_gbif[k]
      res <- count_csv_rows_cached(merge_std$gbif_file[i], cc, cache_path)
      merge_std$n_gbif_after_strict_dedup[i] <- res$n
      cc <- res$cc
      if (k %% 5 == 0) cat("[status_table] GBIF back-fill: ", k, "/", length(needs_gbif), "\n", sep = "")
    }
  }
  
  if (length(needs_nbn) > 0) {
    for (k in seq_along(needs_nbn)) {
      i <- needs_nbn[k]
      res <- count_csv_rows_cached(merge_std$nbn_file[i], cc, cache_path)
      merge_std$n_nbn_after_strict_dedup[i] <- res$n
      cc <- res$cc
      if (k %% 10 == 0) cat("[status_table] NBN back-fill: ", k, "/", length(needs_nbn), "\n", sep = "")
    }
  }
  
  if (length(needs_out) > 0) {
    for (k in seq_along(needs_out)) {
      i <- needs_out[k]
      res <- count_table_rows_cached(merge_std$out_file[i], cc, cache_path)
      merge_std$n_after_merge[i] <- res$n
      cc <- res$cc
      if (k %% 20 == 0) cat("[status_table] merged back-fill: ", k, "/", length(needs_out), "\n", sep = "")
    }
  }
  
  miss_pre <- which(is.na(merge_std$n_total_after_strict_dedup) &
                      !is.na(merge_std$n_gbif_after_strict_dedup) &
                      !is.na(merge_std$n_nbn_after_strict_dedup))
  if (length(miss_pre) > 0) {
    merge_std$n_total_after_strict_dedup[miss_pre] <-
      as.integer(merge_std$n_gbif_after_strict_dedup[miss_pre] + merge_std$n_nbn_after_strict_dedup[miss_pre])
  }
  
  merge_std[is.na(gbif_downloaded) & !is.na(gbif_file) & file.exists(gbif_file), gbif_downloaded := TRUE]
  merge_std[is.na(nbn_downloaded)  & !is.na(nbn_file)  & file.exists(nbn_file),  nbn_downloaded  := TRUE]
}

# ---- Stage 03 QC flagging -----------------------------------------------------
qc_dt <- safe_read_csv_dt(qc_runlog)
qc_latest <- latest_by_slug(qc_dt)

qc_std <- NULL
if (!is.null(qc_latest) && nrow(qc_latest) > 0) {
  c_status <- pick_col(qc_latest, c("status"))
  c_note   <- pick_col(qc_latest, c("note"))
  c_ts     <- pick_col(qc_latest, c("timestamp_utc", "timestamp"))
  c_n_in   <- pick_col(qc_latest, c("n_in", "n", "n_read"))
  
  qc_std <- qc_latest[, .(
    slug = slug,
    n_after_qc_flagging = if (!is.null(c_n_in)) suppressWarnings(as.integer(get(c_n_in))) else NA_integer_,
    stage03_qc_status = if (!is.null(c_status)) as.character(get(c_status)) else NA_character_,
    stage03_qc_time_utc = if (!is.null(c_ts)) as.character(get(c_ts)) else NA_character_,
    stage03_qc_note = if (!is.null(c_note)) as.character(get(c_note)) else NA_character_
  )]
}

# ---- Stage 04 filtering -------------------------------------------------------
filter_dt <- safe_read_csv_dt(filter_runlog)
filter_latest <- latest_by_slug(filter_dt)

filter_std <- NULL
if (!is.null(filter_latest) && nrow(filter_latest) > 0) {
  c_status <- pick_col(filter_latest, c("status"))
  c_note   <- pick_col(filter_latest, c("note"))
  c_ts     <- pick_col(filter_latest, c("timestamp_utc", "timestamp"))
  c_policy <- pick_col(filter_latest, c("policy_id", "policy"))
  c_n_in   <- pick_col(filter_latest, c("n_in", "n_read", "n_before"))
  c_n_out  <- pick_col(filter_latest, c("n_out", "n_after", "n_filtered"))
  
  filter_std <- filter_latest[, .(
    slug = slug,
    filter_policy_id = if (!is.null(c_policy)) as.character(get(c_policy)) else NA_character_,
    n_before_policy_filter = if (!is.null(c_n_in)) suppressWarnings(as.integer(get(c_n_in))) else NA_integer_,
    n_after_policy_filter  = if (!is.null(c_n_out)) suppressWarnings(as.integer(get(c_n_out))) else NA_integer_,
    stage04_filter_status = if (!is.null(c_status)) as.character(get(c_status)) else NA_character_,
    stage04_filter_time_utc = if (!is.null(c_ts)) as.character(get(c_ts)) else NA_character_,
    stage04_filter_note = if (!is.null(c_note)) as.character(get(c_note)) else NA_character_
  )]
  
  filter_std[, kept_pct_after_policy := ifelse(
    !is.na(n_before_policy_filter) & n_before_policy_filter > 0L & !is.na(n_after_policy_filter),
    round(100 * n_after_policy_filter / n_before_policy_filter, 1),
    NA_real_
  )]
  filter_std[, dropped_by_policy := ifelse(
    !is.na(n_before_policy_filter) & !is.na(n_after_policy_filter),
    n_before_policy_filter - n_after_policy_filter,
    NA_integer_
  )]
}

# ---- Stage 05 gridding --------------------------------------------------------
grid_std <- NULL
grid_source_path <- if (!is.na(grid_runlog)) grid_runlog else grid_summary
gs <- safe_read_csv_dt(grid_source_path)

if (!is.null(gs) && nrow(gs) > 0) {
  c_slug <- pick_col(gs, c("slug", "species_slug"))
  if (!is.null(c_slug)) setnames(gs, c_slug, "slug")
  
  c_status <- pick_col(gs, c("status"))
  c_note   <- pick_col(gs, c("note"))
  c_ts     <- pick_col(gs, c("timestamp_utc", "timestamp"))
  
  c_presence_cells <- pick_col(gs, c("n_presence_cells", "presence_cells", "n_cells_presence"))
  c_gee_points     <- pick_col(gs, c("gee_points", "n_gee_points", "n_points_gee"))
  
  policy_tag_guess <- tryCatch(basename(dirname(grid_source_path)), error = function(e) NA_character_)
  
  grid_std <- gs[, .(
    slug = slug,
    stage05_grid_policy_tag = policy_tag_guess,
    n_gridded_presence_cells = if (!is.null(c_presence_cells)) suppressWarnings(as.integer(get(c_presence_cells))) else NA_integer_,
    n_gridded_points_for_gee = if (!is.null(c_gee_points)) suppressWarnings(as.integer(get(c_gee_points))) else NA_integer_,
    stage05_grid_status = if (!is.null(c_status)) as.character(get(c_status)) else NA_character_,
    stage05_grid_time_utc = if (!is.null(c_ts)) as.character(get(c_ts)) else NA_character_,
    stage05_grid_note = if (!is.null(c_note)) as.character(get(c_note)) else NA_character_,
    stage05_grid_source_path = grid_source_path
  )]
}

# ---- Assemble -----------------------------------------------------------------
out <- copy(base)

if (!is.null(merge_std))  out <- merge(out, merge_std,  by = "slug", all.x = TRUE)
if (!is.null(qc_std))     out <- merge(out, qc_std,     by = "slug", all.x = TRUE)
if (!is.null(filter_std)) out <- merge(out, filter_std, by = "slug", all.x = TRUE)
if (!is.null(grid_std))   out <- merge(out, grid_std,   by = "slug", all.x = TRUE)

out[, gbif := yesno(gbif_downloaded)]
out[, nbn  := yesno(nbn_downloaded)]

# Placeholders for later stages / manual updates
out[, gee_assets_exist := ""]
out[, gee_model_run := ""]
out[, gee_model_refinement_status := ""]
out[, records_lost_sensitive_species := ""]

# ---- Column order -------------------------------------------------------------
keep_cols <- c(
  "species_binomial",
  "slug",
  
  "raw_group_dir_used",
  "gbif",
  "n_gbif_after_strict_dedup",
  "nbn",
  "n_nbn_after_strict_dedup",
  
  "n_total_after_strict_dedup",
  "n_dropped_between_sources",
  "n_after_merge",
  "stage02_merge_status",
  
  "n_after_qc_flagging",
  "stage03_qc_status",
  
  "filter_policy_id",
  "n_before_policy_filter",
  "n_after_policy_filter",
  "kept_pct_after_policy",
  "dropped_by_policy",
  "stage04_filter_status",
  
  "stage05_grid_policy_tag",
  "n_gridded_presence_cells",
  "n_gridded_points_for_gee",
  "stage05_grid_status",
  
  "gee_assets_exist",
  "gee_model_run",
  "gee_model_refinement_status",
  "records_lost_sensitive_species",
  
  "stage02_merge_note",
  "stage03_qc_note",
  "stage04_filter_note",
  "stage05_grid_note",
  "stage05_grid_source_path"
)

for (nm in keep_cols) if (!nm %in% names(out)) out[, (nm) := NA]
out <- out[, ..keep_cols]

# ---- Write output -------------------------------------------------------------
fwrite(out, out_path, na = "")

cat("\n============================================================\n")
cat("Species pipeline status table written\n")
cat("============================================================\n")
cat("Output: ", out_path, "\n", sep = "")
cat("Count cache: ", cache_path, "\n", sep = "")
cat("============================================================\n\n")