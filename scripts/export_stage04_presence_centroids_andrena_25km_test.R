# scripts/export_stage04_presence_centroids_andrena_25km_test.R ------------------
# Stage 04 -> Earth Engine helper:
# Export presence-cell centroids for Andrena fulva from the Stage 04 grid output.
#
# Expected inputs:
# - Stage 04 output folder exists, e.g.
#   data/processed/04_grid/grid25km_test_landmask_europe_bbox/andrena_fulva/...
# - A grid table file exists in that folder (CSV/Parquet/RDS) containing:
#   presence, lon_center, lat_center, cell_id (optional), n_points_in_cell (optional)
#
# Output:
# - A CSV with one row per presence cell:
#   lon, lat, cell_id, n_points_in_cell
#
# Checks printed:
# - number of presence cells (expect ~1565 for your current run)
# - lon/lat ranges
# - a few example rows

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Find this script’s directory (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run this via source('.../scripts/export_stage04_presence_centroids_andrena_25km_test.R') from a file."
  )
}

script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

# ==============================================================================
# CONTROL PANEL (should match  Stage 04 run)
# ==============================================================================
out_stage  <- "04_grid"
policy_tag <- "grid25km_test_landmask_europe_bbox"
cell_km    <- 25

species_name <- "Andrena fulva"

# Europe-ish bbox (for sanity checks)
bbox_ll <- list(xmin = -25, xmax = 45, ymin = 34, ymax = 72)

# ==============================================================================
# Helpers
# ==============================================================================
slugify_species <- function(species_name) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(species_name))
  slug <- gsub("^_+|_+$", "", slug)
  slug
}

guess_grid_table_file <- function(species_dir) {
  # Prefer Parquet, then CSV, then RDS
  candidates <- c(
    list.files(species_dir, pattern = "\\.parquet$", full.names = TRUE, ignore.case = TRUE),
    list.files(species_dir, pattern = "\\.csv$",     full.names = TRUE, ignore.case = TRUE),
    list.files(species_dir, pattern = "\\.rds$",     full.names = TRUE, ignore.case = TRUE)
  )
  if (length(candidates) == 0) return(NA_character_)
  
  # Heuristic: prefer files with "grid" in name, else first candidate.
  grid_like <- candidates[grepl("grid", basename(candidates), ignore.case = TRUE)]
  if (length(grid_like) > 0) return(grid_like[1])
  candidates[1]
}

read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    return(fread(path))
  } else if (ext == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Need the 'arrow' package to read Parquet. Install.packages('arrow') then retry.")
    }
    return(as.data.table(arrow::read_parquet(path)))
  } else if (ext == "rds") {
    obj <- readRDS(path)
    return(as.data.table(obj))
  } else {
    stop("Unsupported file type: ", ext)
  }
}

# ==============================================================================
# Locate Stage 04 grid output for Andrena
# ==============================================================================
slug <- slugify_species(species_name)

species_dir <- file.path("data", "processed", out_stage, policy_tag, slug)
if (!dir.exists(species_dir)) {
  stop(
    "Can't find species Stage 04 directory:\n  ", normalizePath(species_dir, winslash = "/"),
    "\nHave you run scripts/grid_stage03_25km_test.R successfully for this policy_tag?"
  )
}

grid_file <- guess_grid_table_file(species_dir)
if (is.na(grid_file) || !file.exists(grid_file)) {
  stop(
    "Can't find a grid table file in:\n  ", normalizePath(species_dir, winslash = "/"),
    "\nExpected at least one .parquet, .csv, or .rds."
  )
}

message("[export] Reading grid table: ", normalizePath(grid_file, winslash = "/"))
g <- read_any_table(grid_file)

# ==============================================================================
# Validate columns (allow small name variations)
# ==============================================================================
# Required
col_presence <- intersect(names(g), c("presence", "pres", "is_presence"))
col_lon      <- intersect(names(g), c("lon_center", "lon_cent", "lon", "longitude"))
col_lat      <- intersect(names(g), c("lat_center", "lat_cent", "lat", "latitude"))

if (length(col_presence) == 0) stop("Can't find a 'presence' column in the grid table.")
if (length(col_lon) == 0)      stop("Can't find a lon centre column (e.g. 'lon_center').")
if (length(col_lat) == 0)      stop("Can't find a lat centre column (e.g. 'lat_center').")

col_presence <- col_presence[1]
col_lon <- col_lon[1]
col_lat <- col_lat[1]

# Optional
col_cell_id <- intersect(names(g), c("cell_id", "grid_id", "id"))
col_npts    <- intersect(names(g), c("n_points_in_cell", "n_points", "count", "n"))

col_cell_id <- if (length(col_cell_id) > 0) col_cell_id[1] else NA_character_
col_npts    <- if (length(col_npts) > 0)    col_npts[1]    else NA_character_

# ==============================================================================
# Export presence cells
# ==============================================================================
# Coerce presence to numeric-ish
g[, (col_presence) := as.integer(get(col_presence))]

pres <- g[get(col_presence) == 1, .(
  lon = as.numeric(get(col_lon)),
  lat = as.numeric(get(col_lat)),
  cell_id = if (!is.na(col_cell_id)) as.character(get(col_cell_id)) else NA_character_,
  n_points_in_cell = if (!is.na(col_npts)) as.integer(get(col_npts)) else NA_integer_
)]

# Drop any missing coords
pres <- pres[is.finite(lon) & is.finite(lat)]

# Write output next to the grid files (easy to find)
out_csv <- file.path(species_dir, sprintf("presence_centroids_%dkm_%s.csv", cell_km, slug))
fwrite(pres, out_csv)

# ==============================================================================
# Checks / printout
# ==============================================================================
message("\n============================================================")
message("Stage 04 -> EE export done")
message("Species: ", species_name, " (", slug, ")")
message("Wrote:   ", normalizePath(out_csv, winslash = "/"))
message("Rows (presence cells): ", nrow(pres))
message(sprintf("Lon range: %.3f .. %.3f", min(pres$lon), max(pres$lon)))
message(sprintf("Lat range: %.3f .. %.3f", min(pres$lat), max(pres$lat)))

# bbox sanity check
oob <- pres[lon < bbox_ll$xmin | lon > bbox_ll$xmax | lat < bbox_ll$ymin | lat > bbox_ll$ymax]
message("Out-of-bbox rows (should be 0 or tiny): ", nrow(oob))

message("\nFirst 10 rows:")
print(pres[1:min(10, .N)])

message("============================================================\n")
