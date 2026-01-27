# scripts/rasterise_stage03_25km_test.R -----------------------------------------
# InfluentialSpecies - Ad hoc rasterisation / gridding from Stage 03 outputs
#
# Purpose:
#   Convert Stage 03 filtered point occurrences into a 25 x 25 km grid for
#   presence/absence modelling (MaxEnt-style background + presences).
#
# Key decisions (current "first-pass" defaults):
#   - Grid CRS: ETRS89 / LAEA Europe (EPSG:3035) for true km-sized cells.
#   - Extent: a simple Europe-ish bounding box (NOT a formal GBIF "continent=EUROPE").
#     We keep this as a pragmatic default and record it here for later revision.
#   - Mask: keep LAND cells only (cell centre falls on land polygon).
#   - Output: a per-species grid table (presence + count). Optional GeoTIFFs if terra is installed.
#
# Inputs:
#   data/processed/03_filtered/<slug>/occ_<slug>__filtered.(parquet|rds)
#
# Outputs (suggested):
#   data/processed/04_grid25km/<policy_tag>/grid25km_land_cells.(parquet|csv)
#   data/processed/04_grid25km/<policy_tag>/<slug>/occ_<slug>__grid25km.(parquet|csv)
#   data/processed/04_grid25km/<policy_tag>/<slug>/presence_25km.tif (optional)
#   data/processed/04_grid25km/<policy_tag>/<slug>/count_25km.tif    (optional)

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Find this script’s directory (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run this via source('.../scripts/rasterise_stage03_25km_test.R') from a file, not by copy/paste."
  )
}

script_dir <- dirname(normalizePath(this_file))
repo_root  <- normalizePath(file.path(script_dir, ".."))
setwd(repo_root)

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# ---- Species list -------------------------------------------------------------
species_names <- c(
  "Myrmica sabuleti",
  "Myrmica scabrinodis",
  "Andrena fulva",
  "Sorex araneus",
  "Leptothorax acervorum",
  "Emberiza schoeniclus"
)

# ---- Where to read Stage 03 from ---------------------------------------------
in_root <- file.path("data", "processed", "03_filtered")

# ---- Where to write gridded outputs ------------------------------------------
# "policy_tag" is just a folder label for organisation (does not need to match Stage 03 policy_id).
policy_tag <- "grid25km_test_landmask_europe_bbox"
out_root <- file.path("data", "processed", "04_grid25km", policy_tag)

# ---- Grid settings ------------------------------------------------------------
cell_km <- 25
cell_m  <- cell_km * 1000

# CRS for a Europe-wide equal-area metric grid
crs_grid <- "EPSG:3035"  # ETRS89 / LAEA Europe

# Europe-ish extent as a bounding box in lon/lat (WGS84).
# NOTE: This is a pragmatic modelling extent for the first pass.
# It is NOT guaranteed to match GBIF's interpreted "continent=EUROPE" logic.
bbox_ll <- list(
  xmin = -25,  # lon
  xmax =  45,
  ymin =  34,  # lat
  ymax =  72
)

# Land mask: keep only cells whose CENTRE falls on land (fast heuristic).
use_land_mask <- TRUE

# Optional: write GeoTIFF rasters if terra is installed
write_geotiff <- TRUE

# ==============================================================================
# HELPERS
# ==============================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

slugify_species <- function(species_name) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(species_name))
  slug <- gsub("^_+|_+$", "", slug)
  slug
}

.ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

.read_stage03_base <- function(repo_root, slug, in_root) {
  base <- file.path(repo_root, in_root, slug, paste0("occ_", slug, "__filtered"))
  p_parq <- paste0(base, ".parquet")
  p_rds  <- paste0(base, ".rds")
  
  if (file.exists(p_parq)) return(list(path = p_parq, fmt = "parquet", base = base))
  if (file.exists(p_rds))  return(list(path = p_rds,  fmt = "rds",     base = base))
  list(path = NA_character_, fmt = NA_character_, base = base)
}

.safe_num <- function(x) suppressWarnings(as.numeric(x))

# ---- Read only what we need (coordinates) ------------------------------------
.read_coords_stage03 <- function(path, fmt) {
  if (is.na(path) || !nzchar(path)) return(NULL)
  
  if (fmt == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Parquet input found but 'arrow' is not installed: install.packages('arrow')")
    }
    
    # Read minimal columns if they exist
    # (Stage 03 typically has lon/lat; fall back to decimalLongitude/decimalLatitude if needed)
    schema_names <- tryCatch(arrow::read_schema(path)$names, error = function(e) character())
    want <- intersect(
      c("lon", "lat", "decimalLongitude", "decimalLatitude", "source"),
      schema_names
    )
    
    dt <- arrow::read_parquet(path, col_select = want)
    setDT(dt)
    return(dt)
  }
  
  if (fmt == "rds") {
    dt <- readRDS(path)
    setDT(dt)
    return(dt)
  }
  
  stop("Unknown format: ", fmt)
}

# ---- Project lon/lat -> EPSG:3035 without building sf geometries (fast) -------
.project_ll_to_grid <- function(lon, lat, crs_to = crs_grid) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Need 'sf' for projection: install.packages('sf')")
  }
  
  # sf_project expects a matrix of (x=lon, y=lat) in EPSG:4326
  pts <- cbind(lon, lat)
  
  if ("sf_project" %in% getNamespaceExports("sf")) {
    xy <- sf::sf_project(from = "EPSG:4326", to = crs_to, pts = pts)
    return(list(x = xy[, 1], y = xy[, 2]))
  }
  
  # Fallback if sf_project is unavailable (older sf): use st_transform (heavier)
  sf_pts <- sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326)
  sf_pts <- sf::st_transform(sf_pts, crs_to)
  xy <- sf::st_coordinates(sf_pts)
  list(x = xy[, 1], y = xy[, 2])
}

# ---- Build a snapped grid extent + indices -----------------------------------
.make_snapped_extent <- function(xmin, xmax, ymin, ymax, cell_m) {
  xmin_s <- floor(xmin / cell_m) * cell_m
  ymin_s <- floor(ymin / cell_m) * cell_m
  xmax_s <- ceiling(xmax / cell_m) * cell_m
  ymax_s <- ceiling(ymax / cell_m) * cell_m
  
  ncol <- as.integer((xmax_s - xmin_s) / cell_m)
  nrow <- as.integer((ymax_s - ymin_s) / cell_m)
  
  list(
    xmin = xmin_s, xmax = xmax_s,
    ymin = ymin_s, ymax = ymax_s,
    ncol = ncol, nrow = nrow
  )
}

# ---- Land polygons (for land masking) ----------------------------------------
.get_land_union <- function(crs_to = crs_grid) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Need 'sf' for land mask: install.packages('sf')")
  }
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
    stop("Need 'rnaturalearth' for land mask: install.packages('rnaturalearth')")
  }
  
  land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  land <- sf::st_make_valid(land)
  land <- sf::st_union(land)  # one (multi)polygon
  land <- sf::st_transform(land, crs_to)
  land
}

# ==============================================================================
# GRID SETUP (Europe bbox -> EPSG:3035 grid -> optional land mask)
# ==============================================================================

.ensure_dir(file.path(repo_root, out_root))

if (!requireNamespace("sf", quietly = TRUE)) {
  stop("This script requires 'sf'. Install it first: install.packages('sf')")
}

# Project bbox corners to grid CRS to get a metric extent
bbox_lon <- c(bbox_ll$xmin, bbox_ll$xmax, bbox_ll$xmax, bbox_ll$xmin)
bbox_lat <- c(bbox_ll$ymin, bbox_ll$ymin, bbox_ll$ymax, bbox_ll$ymax)
bbox_xy  <- .project_ll_to_grid(bbox_lon, bbox_lat, crs_to = crs_grid)

ext <- .make_snapped_extent(
  xmin = min(bbox_xy$x), xmax = max(bbox_xy$x),
  ymin = min(bbox_xy$y), ymax = max(bbox_xy$y),
  cell_m = cell_m
)

# Build a table of ALL grid cell centres (full extent)
grid_dt <- CJ(
  col = 1:ext$ncol,
  row = 1:ext$nrow
)
grid_dt[, `:=`(
  x_center = ext$xmin + (col - 0.5) * cell_m,
  y_center = ext$ymin + (row - 0.5) * cell_m
)]
grid_dt[, cell_id := sprintf("r%05d_c%05d", row, col)]

# Land mask (centre-point heuristic)
if (isTRUE(use_land_mask)) {
  cat("\n[grid] Building land mask (cell-centre heuristic via Natural Earth)...\n")
  land_union <- .get_land_union(crs_to = crs_grid)
  
  centers_sf <- sf::st_as_sf(
    grid_dt[, .(cell_id, x_center, y_center)],
    coords = c("x_center", "y_center"),
    crs = crs_grid
  )
  
  on_land <- as.logical(sf::st_intersects(centers_sf, land_union, sparse = FALSE)[, 1])
  grid_dt[, on_land := on_land]
  grid_land <- grid_dt[on_land == TRUE]
} else {
  grid_dt[, on_land := TRUE]
  grid_land <- grid_dt
}

cat("\n============================================================\n")
cat("Grid setup\n")
cat("CRS:", crs_grid, "\n")
cat("Cell size:", cell_km, "km\n")
cat("Extent (snapped, metres):\n")
cat("  xmin:", ext$xmin, " xmax:", ext$xmax, "\n")
cat("  ymin:", ext$ymin, " ymax:", ext$ymax, "\n")
cat("Cells (all):", nrow(grid_dt), " | land:", nrow(grid_land), "\n")
cat("============================================================\n")

# Save the land grid once (useful for debugging and later reuse)
grid_out_base <- file.path(repo_root, out_root, "grid25km_land_cells")
if (requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(as.data.frame(grid_land), paste0(grid_out_base, ".parquet"))
} else {
  fwrite(grid_land, paste0(grid_out_base, ".csv"))
}

# ==============================================================================
# PER-SPECIES GRIDDING
# ==============================================================================

# Optional terra rasters
terra_ok <- isTRUE(write_geotiff) && requireNamespace("terra", quietly = TRUE)
if (isTRUE(write_geotiff) && !terra_ok) {
  message("[raster] 'terra' not installed; skipping GeoTIFF outputs (tables will still be written).")
}

.assign_to_cells <- function(x, y, ext, cell_m) {
  # returns row/col (1-based) and a logical in_bounds
  col <- as.integer(floor((x - ext$xmin) / cell_m) + 1L)
  row <- as.integer(floor((y - ext$ymin) / cell_m) + 1L)
  
  inb <- !is.na(col) & !is.na(row) &
    col >= 1L & col <= ext$ncol &
    row >= 1L & row <= ext$nrow
  
  list(row = row, col = col, in_bounds = inb)
}

for (sp in species_names) {
  slug <- slugify_species(sp)
  info <- .read_stage03_base(repo_root, slug, in_root)
  
  cat("\n------------------------------------------------------------\n")
  cat("Species:", sp, "\n")
  cat("Slug:", slug, "\n")
  
  if (is.na(info$path)) {
    cat("[skip] No Stage 03 file found.\n")
    next
  }
  
  cat("Reading Stage 03:", info$path, "\n")
  dt <- .read_coords_stage03(info$path, info$fmt)
  if (is.null(dt) || nrow(dt) == 0) {
    cat("[skip] Empty dataset.\n")
    next
  }
  
  # Standardise coordinate columns
  if (!("lon" %in% names(dt)) || !("lat" %in% names(dt))) {
    if (all(c("decimalLongitude", "decimalLatitude") %in% names(dt))) {
      dt[, `:=`(
        lon = .safe_num(decimalLongitude),
        lat = .safe_num(decimalLatitude)
      )]
    }
  } else {
    dt[, `:=`(
      lon = .safe_num(lon),
      lat = .safe_num(lat)
    )]
  }
  
  dt <- dt[!is.na(lon) & !is.na(lat)]
  if (nrow(dt) == 0) {
    cat("[skip] No usable coordinates.\n")
    next
  }
  
  # Project to grid CRS (fast path)
  xy <- .project_ll_to_grid(dt$lon, dt$lat, crs_to = crs_grid)
  
  # Assign points to cells
  idx <- .assign_to_cells(xy$x, xy$y, ext = ext, cell_m = cell_m)
  dt[, `:=`(row = idx$row, col = idx$col, in_bounds = idx$in_bounds)]
  dt <- dt[in_bounds == TRUE]
  
  if (nrow(dt) == 0) {
    cat("[warn] All points fell outside the grid extent (check bbox_ll).\n")
    next
  }
  
  # Count points per cell
  counts <- dt[, .(n_points_in_cell = .N), by = .(row, col)]
  counts[, cell_id := sprintf("r%05d_c%05d", row, col)]
  
  # Join onto land grid to get a full land-background table
  out <- merge(
    grid_land[, .(cell_id, row, col, x_center, y_center, on_land)],
    counts[, .(cell_id, n_points_in_cell)],
    by = "cell_id",
    all.x = TRUE
  )
  out[is.na(n_points_in_cell), n_points_in_cell := 0L]
  out[, presence := as.integer(n_points_in_cell > 0L)]
  
  # Optional lon/lat centre (handy for non-projected tools)
  # Convert centres back to lon/lat (WGS84)
  back <- sf::sf_project(from = crs_grid, to = "EPSG:4326", pts = as.matrix(out[, .(x_center, y_center)]))
  out[, `:=`(
    lon_center = back[, 1],
    lat_center = back[, 2]
  )]
  
  # Write per-species table
  sp_dir <- file.path(repo_root, out_root, slug)
  .ensure_dir(sp_dir)
  
  out_base <- file.path(sp_dir, paste0("occ_", slug, "__grid25km"))
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(as.data.frame(out), paste0(out_base, ".parquet"))
  } else {
    fwrite(out, paste0(out_base, ".csv"))
  }
  
  # Optional GeoTIFFs (presence + counts)
  if (terra_ok) {
    r <- terra::rast(
      ncols = ext$ncol, nrows = ext$nrow,
      xmin = ext$xmin, xmax = ext$xmax,
      ymin = ext$ymin, ymax = ext$ymax,
      crs  = crs_grid
    )
    
    # Build full-grid vectors (NA outside land, values on land)
    # We'll start NA everywhere, then fill land cells using row/col -> cell index
    ncell <- terra::ncell(r)
    v_presence <- rep(NA_integer_, ncell)
    v_count    <- rep(NA_integer_, ncell)
    
    # Convert (row,col) to terra cell index:
    # terra uses row=1 at top; our "row" counts from ymin upwards.
    # So convert: terra_row = (nrow - row + 1)
    terra_row <- ext$nrow - out$row + 1L
    cell_idx  <- terra::cellFromRowCol(r, terra_row, out$col)
    
    # Fill land cells
    v_presence[cell_idx] <- out$presence
    v_count[cell_idx]    <- as.integer(out$n_points_in_cell)
    
    terra::values(r) <- v_presence
    terra::writeRaster(r, file.path(sp_dir, "presence_25km.tif"), overwrite = TRUE)
    
    terra::values(r) <- v_count
    terra::writeRaster(r, file.path(sp_dir, "count_25km.tif"), overwrite = TRUE)
  }
  
  cat("[done] Land cells:", nrow(out), " | presence cells:", sum(out$presence), " | points total:", sum(out$n_points_in_cell), "\n")
}

cat("\n============================================================\n")
cat("Done.\n")
cat("Outputs:", file.path(repo_root, out_root), "\n")
cat("Reminder: extent is a simple Europe-ish bbox (see bbox_ll). Revisit later if needed.\n")
cat("============================================================\n")
