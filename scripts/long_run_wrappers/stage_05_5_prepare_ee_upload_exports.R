#!/usr/bin/env Rscript

# ==============================================================================
# InfluentialSpecies — Stage 05.5 — Prepare Stage 05 presence-point exports
# for manual Earth Engine upload
# ==============================================================================
#
# Purpose
#   Read the per-species Stage 05 presence-point CSVs, standardise them into one
#   Earth Engine-ready schema, split them into upload shards, and write:
#
#   1) a per-species summary table
#   2) a per-shard summary table
#   3) shard CSVs ready for manual upload into Earth Engine
#
# Why this stage exists
#   Stage 05 produces one 1 km representative-point CSV per species. That is
#   convenient for gridding and checking individual species, but awkward for
#   Earth Engine ingestion across the full species set. This stage creates one
#   standard export package with consistent columns and manageable batch sizes.
#
# Earth Engine upload shape used here
#   - CSV tables
#   - point geometry inferred from numeric longitude / latitude columns
#   - EPSG:4326 coordinates
#   - multiple shards with identical schema
#
# Main outputs
#   data/processed/05_5_ee_export/<export_tag>/
#     _summary_species.csv
#     _summary_shards.csv
#     shards/
#
# Notes
#   - This script prepares the final R-side handoff only.
#   - Upload happens manually in the Earth Engine Assets tab.
#   - The export uses column names "longitude" and "latitude" so Earth Engine
#     can detect the point columns without manual remapping.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# Helpers
# ==============================================================================

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x

slugify_species <- function(species_name) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(species_name))
  slug <- gsub("^_+|_+$", "", slug)
  slug
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))
safe_int <- function(x) suppressWarnings(as.integer(x))
safe_chr <- function(x) as.character(x)

find_repo_root <- function(start_dir) {
  markers <- c(".git", "R", "data", "README.md", "scripts")
  d <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  for (i in 1:25) {
    if (sum(file.exists(file.path(d, markers))) >= 3L) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

detect_coord_pair <- function(nms) {
  pairs <- list(
    c("longitude", "latitude"),
    c("lon", "lat"),
    c("decimalLongitude", "decimalLatitude"),
    c("decimal_longitude", "decimal_latitude"),
    c("x", "y"),
    c("X", "Y")
  )
  for (p in pairs) {
    if (all(p %in% nms)) return(p)
  }
  NULL
}

first_present <- function(dt, candidates, default = NA) {
  hit <- intersect(candidates, names(dt))
  if (length(hit) == 0L) return(rep(default, nrow(dt)))
  dt[[hit[1]]]
}

assert_required_dir <- function(path, what) {
  if (!dir.exists(path)) stop(what, " not found: ", path)
}

assert_required_file <- function(path, what) {
  if (!file.exists(path)) stop(what, " not found: ", path)
}

# ==============================================================================
# Resolve script and repo paths
# ==============================================================================

script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_path) || !nzchar(script_path)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run via source('.../stage_05_5_prepare_ee_upload_exports.R') or Rscript."
  )
}

script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
repo_root  <- find_repo_root(script_dir)

# ==============================================================================
# Control panel
# ==============================================================================

# ---- Canonical handoff from Stage 05 -----------------------------------------
policy_tag <- "grid1km_first_observed_landmask_europe_bbox"
stage05_root <- file.path(repo_root, "data", "processed", "05_grid", policy_tag)

# ---- Output identity ----------------------------------------------------------
out_stage  <- "05_5_ee_export"
export_tag <- paste0(policy_tag, "__ee_manual_upload")
out_root   <- file.path(repo_root, "data", "processed", out_stage, export_tag)

shards_dir <- file.path(out_root, "shards")

# ---- Species list -------------------------------------------------------------
species_list_path <- file.path(repo_root, "data", "_meta", "species_list_binomial.csv")

# ---- Export structure ---------------------------------------------------------
# One row = one Stage 05 representative point for one occupied 1 km cell.
# Shards are sized by row count to keep manual uploads manageable.
max_rows_per_shard <- 1000000L

# ---- Output writing -----------------------------------------------------------
overwrite <- TRUE

# ---- Run behaviour ------------------------------------------------------------
strict_missing_species <- TRUE
verbose <- TRUE

# ==============================================================================
# Pre-flight
# ==============================================================================

assert_required_dir(stage05_root, "Stage 05 policy folder")
assert_required_file(species_list_path, "Canonical species list")

ensure_dir(out_root)
ensure_dir(shards_dir)

species_names <- readLines(species_list_path, warn = FALSE, encoding = "UTF-8")
species_names <- species_names[!is.na(species_names)]
species_names <- trimws(gsub("\\s+", " ", species_names))
species_names <- species_names[nzchar(species_names)]
species_names <- species_names[tolower(species_names) != "binomial"]

if (length(species_names) != 100L) {
  stop(
    "Canonical species list did not read as 100 names (got ", length(species_names), "). ",
    "Check: ", species_list_path
  )
}

cat("\n============================================================\n")
cat("Stage 05.5 — Earth Engine manual-upload export\n")
cat("============================================================\n")
cat("repo_root:           ", repo_root, "\n", sep = "")
cat("stage05_root:        ", stage05_root, "\n", sep = "")
cat("out_root:            ", out_root, "\n", sep = "")
cat("species_n:           ", length(species_names), "\n", sep = "")
cat("policy_tag:          ", policy_tag, "\n", sep = "")
cat("max_rows_per_shard:  ", format(max_rows_per_shard, big.mark = ","), "\n", sep = "")
cat("============================================================\n\n")

# ==============================================================================
# Read and standardise per-species Stage 05 presence files
# ==============================================================================

species_chunks <- vector("list", length(species_names))
species_summary <- vector("list", length(species_names))
missing_species <- character()

for (i in seq_along(species_names)) {
  sp   <- species_names[i]
  slug <- slugify_species(sp)
  
  in_csv <- file.path(stage05_root, slug, paste0("presence_points_1km_", slug, ".csv"))
  
  if (isTRUE(verbose)) {
    cat("[", sprintf("%03d", i), "/", length(species_names), "] ", sp, " -> ", slug, "\n", sep = "")
  }
  
  if (!file.exists(in_csv)) {
    missing_species <- c(missing_species, slug)
    
    species_summary[[i]] <- data.table(
      species = sp,
      slug = slug,
      input_csv = in_csv,
      found = FALSE,
      n_rows = 0L,
      n_bad_coords = 0L,
      min_longitude = NA_real_,
      max_longitude = NA_real_,
      min_latitude = NA_real_,
      max_latitude = NA_real_,
      status = "missing",
      note = "Stage 05 presence_points CSV not found"
    )
    
    if (isTRUE(verbose)) cat("  [missing] ", in_csv, "\n", sep = "")
    next
  }
  
  dt <- fread(in_csv, showProgress = FALSE, encoding = "UTF-8")
  
  if (nrow(dt) == 0L) {
    species_chunks[[i]] <- data.table()
    
    species_summary[[i]] <- data.table(
      species = sp,
      slug = slug,
      input_csv = in_csv,
      found = TRUE,
      n_rows = 0L,
      n_bad_coords = 0L,
      min_longitude = NA_real_,
      max_longitude = NA_real_,
      min_latitude = NA_real_,
      max_latitude = NA_real_,
      status = "empty",
      note = "Stage 05 presence_points CSV exists but contains zero rows"
    )
    
    if (isTRUE(verbose)) cat("  [empty]\n")
    next
  }
  
  coord_pair <- detect_coord_pair(names(dt))
  if (is.null(coord_pair)) {
    stop(
      "Couldn't find coordinate columns in: ", in_csv, "\n",
      "Available columns: ", paste(names(dt), collapse = ", ")
    )
  }
  
  x_in <- coord_pair[1]
  y_in <- coord_pair[2]
  
  dt[, `:=`(
    longitude = safe_num(get(x_in)),
    latitude  = safe_num(get(y_in))
  )]
  
  bad_coord <- is.na(dt$longitude) | is.na(dt$latitude) |
    dt$longitude < -180 | dt$longitude > 180 |
    dt$latitude  <  -90 | dt$latitude  >  90
  
  n_bad <- as.integer(sum(bad_coord))
  
  if (n_bad > 0L) {
    dt <- dt[!bad_coord]
  }
  
  dt[, `:=`(
    species = sp,
    slug = slug,
    cell_id = safe_chr(first_present(dt, c("cell_id"), default = NA_character_)),
    row = safe_int(first_present(dt, c("row"), default = NA_integer_)),
    col = safe_int(first_present(dt, c("col"), default = NA_integer_)),
    n_points_in_cell = safe_int(first_present(dt, c("n_points_in_cell"), default = NA_integer_)),
    cell_km = safe_int(first_present(dt, c("cell_km"), default = 1L)),
    policy_tag = safe_chr(first_present(dt, c("policy_tag"), default = policy_tag)),
    source_csv = basename(in_csv)
  )]
  
  dt <- dt[, .(
    species,
    slug,
    longitude,
    latitude,
    cell_id,
    row,
    col,
    n_points_in_cell,
    cell_km,
    policy_tag,
    source_csv
  )]
  
  setorder(dt, slug, cell_id, longitude, latitude)
  
  species_chunks[[i]] <- dt
  
  species_summary[[i]] <- data.table(
    species = sp,
    slug = slug,
    input_csv = in_csv,
    found = TRUE,
    n_rows = as.integer(nrow(dt)),
    n_bad_coords = n_bad,
    min_longitude = if (nrow(dt) > 0L) min(dt$longitude, na.rm = TRUE) else NA_real_,
    max_longitude = if (nrow(dt) > 0L) max(dt$longitude, na.rm = TRUE) else NA_real_,
    min_latitude = if (nrow(dt) > 0L) min(dt$latitude, na.rm = TRUE) else NA_real_,
    max_latitude = if (nrow(dt) > 0L) max(dt$latitude, na.rm = TRUE) else NA_real_,
    status = "ok",
    note = if (n_bad > 0L) "Dropped rows with invalid or missing coordinates" else ""
  )
  
  if (isTRUE(verbose)) {
    cat(
      "  [ok] rows=", format(nrow(dt), big.mark = ","),
      if (n_bad > 0L) paste0(" | dropped_bad_coords=", format(n_bad, big.mark = ",")) else "",
      "\n",
      sep = ""
    )
  }
}

summary_species_dt <- rbindlist(species_summary, fill = TRUE, use.names = TRUE)
fwrite(summary_species_dt, file.path(out_root, "_summary_species.csv"))

if (length(missing_species) > 0L && isTRUE(strict_missing_species)) {
  stop(
    "Missing Stage 05 presence files for ", length(missing_species), " species.\n",
    "See: ", file.path(out_root, "_summary_species.csv")
  )
}

all_dt <- rbindlist(species_chunks, fill = TRUE, use.names = TRUE)

if (nrow(all_dt) == 0L) {
  stop("No Stage 05 presence rows were available after standardisation.")
}

setorder(all_dt, slug, cell_id, longitude, latitude)

all_dt[, ee_row_id := sprintf("%09d", seq_len(.N))]
all_dt[, species_point_id := sprintf("%s_%07d", slug, seq_len(.N)), by = slug]

setcolorder(all_dt, c(
  "ee_row_id",
  "species_point_id",
  "species",
  "slug",
  "longitude",
  "latitude",
  "cell_id",
  "row",
  "col",
  "n_points_in_cell",
  "cell_km",
  "policy_tag",
  "source_csv"
))

# ==============================================================================
# Shard assignment
# ==============================================================================

all_dt[, shard_index := as.integer((seq_len(.N) - 1L) %/% max_rows_per_shard) + 1L]
all_dt[, shard_name := sprintf("stage05_presence_shard_%02d", shard_index)]

shard_summary_dt <- all_dt[, .(
  n_rows = .N,
  n_species = uniqueN(slug),
  min_longitude = min(longitude, na.rm = TRUE),
  max_longitude = max(longitude, na.rm = TRUE),
  min_latitude = min(latitude, na.rm = TRUE),
  max_latitude = max(latitude, na.rm = TRUE)
), by = .(shard_index, shard_name)][order(shard_index)]

fwrite(shard_summary_dt, file.path(out_root, "_summary_shards.csv"))

# ==============================================================================
# Write shard CSVs
# ==============================================================================

for (i in seq_len(nrow(shard_summary_dt))) {
  shard_i  <- shard_summary_dt$shard_index[i]
  shard_nm <- shard_summary_dt$shard_name[i]
  shard_csv <- file.path(shards_dir, paste0(shard_nm, ".csv"))
  
  shard_dt <- all_dt[shard_index == shard_i, !c("shard_index", "shard_name")]
  
  if (!file.exists(shard_csv) || isTRUE(overwrite)) {
    fwrite(shard_dt, shard_csv)
  }
}

# ==============================================================================
# Final run summary
# ==============================================================================

species_ok_n      <- summary_species_dt[status == "ok", .N]
species_empty_n   <- summary_species_dt[status == "empty", .N]
species_missing_n <- summary_species_dt[status == "missing", .N]

cat("\n============================================================\n")
cat("Stage 05.5 complete\n")
cat("============================================================\n")
cat("Total rows prepared:      ", format(nrow(all_dt), big.mark = ","), "\n", sep = "")
cat("Species with rows:        ", species_ok_n, "\n", sep = "")
cat("Species empty:            ", species_empty_n, "\n", sep = "")
cat("Species missing:          ", species_missing_n, "\n", sep = "")
cat("Shard count:              ", nrow(shard_summary_dt), "\n", sep = "")
cat("Out root:                 ", out_root, "\n", sep = "")
cat("Species summary:          ", file.path(out_root, "_summary_species.csv"), "\n", sep = "")
cat("Shard summary:            ", file.path(out_root, "_summary_shards.csv"), "\n", sep = "")
cat("Shard CSV folder:         ", shards_dir, "\n", sep = "")
cat("============================================================\n\n")