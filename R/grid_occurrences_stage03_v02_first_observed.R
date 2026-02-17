# R/grid_occurrences_stage03_v02_first_observed.R
#
# ==============================================================================
# Stage 03 filtered occurrences -> regular grid outputs
# ==============================================================================
#
# Key design choice
#   We DO NOT use cell centroids as representative presence locations.
#   Centroids are invented locations and can create a misleading lattice when
#   exported to Google Earth Engine (GEE).
#
# Representative point rule (deterministic)
#   For each grid cell with >=1 observation:
#     - choose ONE observed lon/lat from the underlying points
#     - rule: "first_observed" after sorting by lon then lat within the cell
#
# Outputs per species
#   1) Full land-grid table (parquet preferred; csv fallback)
#      Includes presence (0/1), n_points_in_cell, and representative lon/lat for
#      presence cells (NA for absence cells).
#
#   2) GEE-ready presence-points CSV (ALWAYS written)
#      Presence-only rows with representative observed lon/lat for upload as a
#      point FeatureCollection in GEE.
#      File name: presence_points_<cell_km>km_<slug>.csv
#
# Notes
#   - This engine reads Stage 03 outputs (parquet or rds) and only requires
#     coordinates.
#   - Land mask is cell-centre-on-land (Natural Earth union). Keep stable for now.
#
# Main entry point:
#   grid_stage03_to_grid()
#
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# Helpers
# ==============================================================================

slugify_species <- function(species_name) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(species_name))
  slug <- gsub("^_+|_+$", "", slug)
  slug
}

.ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

.safe_num <- function(x) suppressWarnings(as.numeric(x))

# ---- Read Stage 03 per-species file base -------------------------------------
.read_stage03_base <- function(repo_root, slug, in_root) {
  base <- file.path(repo_root, in_root, slug, paste0("occ_", slug, "__filtered"))
  p_parq <- paste0(base, ".parquet")
  p_rds  <- paste0(base, ".rds")
  
  if (file.exists(p_parq)) return(list(path = p_parq, fmt = "parquet", base = base))
  if (file.exists(p_rds))  return(list(path = p_rds,  fmt = "rds",     base = base))
  list(path = NA_character_, fmt = NA_character_, base = base)
}

# ---- Read only what we need (coordinates) ------------------------------------
.read_coords_stage03 <- function(path, fmt) {
  if (is.na(path) || !nzchar(path)) return(NULL)
  
  if (fmt == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("Parquet input found but 'arrow' is not installed: install.packages('arrow')")
    }
    
    coord_pairs <- list(
      c("lon", "lat"),
      c("longitude", "latitude"),
      c("decimalLongitude", "decimalLatitude"),
      c("decimal_longitude", "decimal_latitude"),
      c("lon_r", "lat_r"),
      c("x", "y"),
      c("X", "Y")
    )
    
    for (p in coord_pairs) {
      dt <- tryCatch({
        if (requireNamespace("tidyselect", quietly = TRUE)) {
          arrow::read_parquet(path, col_select = tidyselect::all_of(p))
        } else {
          arrow::read_parquet(path, col_select = p)
        }
      }, error = function(e) NULL)
      
      if (!is.null(dt)) {
        setDT(dt)
        if (p[1] != "lon") setnames(dt, p[1], "lon")
        if (p[2] != "lat") setnames(dt, p[2], "lat")
        return(dt)
      }
    }
    
    schema_names <- tryCatch(arrow::read_schema(path)$names, error = function(e) character())
    stop(
      "Couldn't find coordinate columns in Stage 03 parquet.\n",
      "Available columns: ", paste(schema_names, collapse = ", ")
    )
  }
  
  if (fmt == "rds") {
    obj <- readRDS(path)
    
    if (inherits(obj, "data.frame")) {
      dt <- as.data.table(obj)
    } else if (inherits(obj, "data.table")) {
      dt <- copy(obj)
    } else {
      stop("Stage 03 RDS is not a data.frame/data.table: ", path)
    }
    
    coord_pairs <- list(
      c("lon", "lat"),
      c("longitude", "latitude"),
      c("decimalLongitude", "decimalLatitude"),
      c("decimal_longitude", "decimal_latitude"),
      c("lon_r", "lat_r"),
      c("x", "y"),
      c("X", "Y")
    )
    
    for (p in coord_pairs) {
      if (all(p %in% names(dt))) {
        dt <- dt[, ..p]
        if (p[1] != "lon") setnames(dt, p[1], "lon")
        if (p[2] != "lat") setnames(dt, p[2], "lat")
        return(dt)
      }
    }
    
    stop(
      "Couldn't find coordinate columns in Stage 03 RDS: ", path, "\n",
      "Available columns: ", paste(names(dt), collapse = ", ")
    )
  }
  
  stop("Unknown input format: ", fmt)
}

# ---- Projection helpers -------------------------------------------------------
.project_ll_to_grid <- function(lon, lat, crs_to) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("This stage requires 'sf'. Install it first: install.packages('sf')")
  }
  
  pts <- sf::st_as_sf(
    data.frame(lon = lon, lat = lat),
    coords = c("lon", "lat"),
    crs = "EPSG:4326"
  )
  pts2 <- sf::st_transform(pts, crs_to)
  xy <- sf::st_coordinates(pts2)
  list(x = xy[, 1], y = xy[, 2])
}

# ---- Extent creation (snapped to cell size) -----------------------------------
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

# ---- Assign points to row/col indices ----------------------------------------
.assign_to_cells <- function(x, y, ext, cell_m) {
  col <- as.integer(floor((x - ext$xmin) / cell_m) + 1L)
  row <- as.integer(floor((y - ext$ymin) / cell_m) + 1L)
  
  in_bounds <- !is.na(row) & !is.na(col) &
    row >= 1L & row <= ext$nrow &
    col >= 1L & col <= ext$ncol
  
  list(row = row, col = col, in_bounds = in_bounds)
}

# ---- Natural Earth land union (fast land-mask heuristic) ----------------------
.get_land_union <- function(crs_to) {
  if (!requireNamespace("rnaturalearth", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE)) {
    stop("Land mask requires rnaturalearth + sf. Install packages 'rnaturalearth' and 'sf'.")
  }
  
  land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  land2 <- sf::st_transform(land, crs_to)
  sf::st_union(land2)
}

# ---- Optional bbox context map (lon/lat) --------------------------------------
plot_europe_bbox_map <- function(bbox_ll, out_file) {
  if (!requireNamespace("rnaturalearth", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE)) {
    stop("Map output requires rnaturalearth + sf.")
  }
  
  land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  bb <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = bbox_ll$xmin, ymin = bbox_ll$ymin, xmax = bbox_ll$xmax, ymax = bbox_ll$ymax),
    crs = sf::st_crs(4326)
  ))
  
  png(out_file, width = 1200, height = 800)
  on.exit(dev.off(), add = TRUE)
  plot(sf::st_geometry(land), col = "grey90", border = "grey60", main = "Europe-ish bbox (lon/lat)")
  plot(bb, add = TRUE, border = "red", lwd = 2)
}

# ==============================================================================
# Main entry point
# ==============================================================================

grid_stage03_to_grid <- function(species_names,
                                 in_root,
                                 out_stage = "04_grid",
                                 policy_tag = "grid_first_observed",
                                 cell_km = 25,
                                 crs_grid = "EPSG:3035",
                                 bbox_ll = list(xmin = -25, xmax = 45, ymin = 34, ymax = 72),
                                 use_land_mask = TRUE,
                                 write_geotiff = TRUE,
                                 plot_bbox_map = FALSE,
                                 bbox_map_filename = "bbox_context_map.png",
                                 verbose = TRUE) {
  
  repo_root <- normalizePath(getwd())
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("This stage requires 'sf'. Install it first: install.packages('sf')")
  }
  
  cell_m <- cell_km * 1000
  
  # Where Stage 04 outputs go
  out_root <- file.path("data", "processed", out_stage, policy_tag)
  .ensure_dir(file.path(repo_root, out_root))
  
  # Optional bbox context map (lon/lat)
  if (isTRUE(plot_bbox_map)) {
    bbox_out <- file.path(repo_root, out_root, bbox_map_filename)
    if (isTRUE(verbose)) cat("[map] Writing bbox context map:", bbox_out, "\n")
    plot_europe_bbox_map(bbox_ll = bbox_ll, out_file = bbox_out)
  }
  
  # ---------------------------------------------------------------------------
  # Grid setup (bbox corners -> EPSG:3035 -> snapped extent)
  # ---------------------------------------------------------------------------
  bbox_lon <- c(bbox_ll$xmin, bbox_ll$xmax, bbox_ll$xmax, bbox_ll$xmin)
  bbox_lat <- c(bbox_ll$ymin, bbox_ll$ymin, bbox_ll$ymax, bbox_ll$ymax)
  bbox_xy  <- .project_ll_to_grid(bbox_lon, bbox_lat, crs_to = crs_grid)
  
  ext <- .make_snapped_extent(
    xmin = min(bbox_xy$x), xmax = max(bbox_xy$x),
    ymin = min(bbox_xy$y), ymax = max(bbox_xy$y),
    cell_m = cell_m
  )
  
  grid_dt <- CJ(col = 1:ext$ncol, row = 1:ext$nrow)
  grid_dt[, `:=`(
    x_center = ext$xmin + (col - 0.5) * cell_m,
    y_center = ext$ymin + (row - 0.5) * cell_m
  )]
  grid_dt[, cell_id := sprintf("r%05d_c%05d", row, col)]
  
  # Land mask (cell-centre heuristic)
  if (isTRUE(use_land_mask)) {
    if (isTRUE(verbose)) cat("[grid] Building land mask (cell-centre heuristic via Natural Earth)...\n")
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
  
  if (isTRUE(verbose)) {
    cat("\n============================================================\n")
    cat("Grid setup\n")
    cat("CRS:", crs_grid, "\n")
    cat("Cell size:", cell_km, "km\n")
    cat("Extent (snapped, metres):\n")
    cat("  xmin:", ext$xmin, " xmax:", ext$xmax, "\n")
    cat("  ymin:", ext$ymin, " ymax:", ext$ymax, "\n")
    cat("Cells (all):", nrow(grid_dt), " | land:", nrow(grid_land), "\n")
    cat("Representative point rule: first_observed (sorted lon/lat within cell)\n")
    cat("GEE output: presence-only points CSV per species\n")
    cat("============================================================\n")
  }
  
  # Save the land grid once (useful for debugging and later reuse)
  grid_out_base <- file.path(repo_root, out_root, paste0("grid", cell_km, "km_land_cells"))
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(as.data.frame(grid_land), paste0(grid_out_base, ".parquet"))
  } else {
    fwrite(grid_land, paste0(grid_out_base, ".csv"))
  }
  
  # Optional rasters for quick GIS / sanity checks
  terra_ok <- isTRUE(write_geotiff) && requireNamespace("terra", quietly = TRUE)
  if (isTRUE(write_geotiff) && !terra_ok && isTRUE(verbose)) {
    message("[raster] 'terra' not installed; skipping GeoTIFF outputs (tables will still be written).")
  }
  
  summary_dt <- data.table(
    species = character(),
    slug = character(),
    stage03_file = character(),
    n_read = integer(),
    n_after_non_na = integer(),
    n_in_bounds = integer(),
    n_points_total = integer(),
    n_presence_cells = integer(),
    gee_points = integer()
  )
  
  for (sp in species_names) {
    
    slug <- slugify_species(sp)
    info <- .read_stage03_base(repo_root, slug, in_root)
    
    if (isTRUE(verbose)) {
      cat("\n------------------------------------------------------------\n")
      cat("Species:", sp, "\n")
      cat("Slug:", slug, "\n")
    }
    
    if (is.na(info$path)) {
      if (isTRUE(verbose)) cat("[skip] No Stage 03 file found.\n")
      summary_dt <- rbind(
        summary_dt,
        data.table(
          species = sp, slug = slug, stage03_file = NA_character_,
          n_read = 0L, n_after_non_na = 0L, n_in_bounds = 0L,
          n_points_total = 0L, n_presence_cells = 0L, gee_points = 0L
        )
      )
      next
    }
    
    if (isTRUE(verbose)) cat("Reading Stage 03:", info$path, "\n")
    dt <- .read_coords_stage03(info$path, info$fmt)
    if (is.null(dt) || nrow(dt) == 0) {
      if (isTRUE(verbose)) cat("[skip] Empty dataset.\n")
      summary_dt <- rbind(
        summary_dt,
        data.table(
          species = sp, slug = slug, stage03_file = info$path,
          n_read = 0L, n_after_non_na = 0L, n_in_bounds = 0L,
          n_points_total = 0L, n_presence_cells = 0L, gee_points = 0L
        )
      )
      next
    }
    
    n_read <- nrow(dt)
    if (isTRUE(verbose)) cat("[read] rows:", n_read, " cols:", paste(names(dt), collapse = ", "), "\n")
    
    dt[, `:=`(
      lon = .safe_num(lon),
      lat = .safe_num(lat)
    )]
    
    dt <- dt[!is.na(lon) & !is.na(lat)]
    n_after_non_na <- nrow(dt)
    if (n_after_non_na == 0) {
      if (isTRUE(verbose)) cat("[skip] No usable coordinates.\n")
      summary_dt <- rbind(
        summary_dt,
        data.table(
          species = sp, slug = slug, stage03_file = info$path,
          n_read = n_read, n_after_non_na = 0L, n_in_bounds = 0L,
          n_points_total = 0L, n_presence_cells = 0L, gee_points = 0L
        )
      )
      next
    }
    
    # Project lon/lat into EPSG:3035 to assign to regular metric grid
    xy <- .project_ll_to_grid(dt$lon, dt$lat, crs_to = crs_grid)
    
    idx <- .assign_to_cells(xy$x, xy$y, ext = ext, cell_m = cell_m)
    dt[, `:=`(row = idx$row, col = idx$col, in_bounds = idx$in_bounds)]
    dt <- dt[in_bounds == TRUE]
    n_in_bounds <- nrow(dt)
    
    if (isTRUE(verbose)) {
      cat("[qc] after non-NA coords:", n_after_non_na, " -> after in-bounds:", n_in_bounds, "\n")
    }
    
    if (n_in_bounds == 0) {
      if (isTRUE(verbose)) cat("[warn] All points fell outside the grid extent (check bbox_ll).\n")
      summary_dt <- rbind(
        summary_dt,
        data.table(
          species = sp, slug = slug, stage03_file = info$path,
          n_read = n_read, n_after_non_na = n_after_non_na, n_in_bounds = 0L,
          n_points_total = 0L, n_presence_cells = 0L, gee_points = 0L
        )
      )
      next
    }
    
    # Count observations per cell
    counts <- dt[, .(n_points_in_cell = .N), by = .(row, col)]
    counts[, cell_id := sprintf("r%05d_c%05d", row, col)]
    
    # Representative observed lon/lat per cell (deterministic)
    # Rule: first row after sorting by lon, then lat, within each (row,col)
    setorder(dt, row, col, lon, lat)
    rep_pt <- dt[, .SD[1], by = .(row, col), .SDcols = c("lon", "lat")]
    rep_pt[, cell_id := sprintf("r%05d_c%05d", row, col)]
    rep_pt <- rep_pt[, .(cell_id, lon, lat)]
    
    # Join onto the land grid:
    # - all land cells present (background space)
    # - counts/presence for this species
    # - representative lon/lat only for presence cells
    out <- merge(
      grid_land[, .(cell_id, row, col, x_center, y_center, on_land)],
      counts[, .(cell_id, n_points_in_cell)],
      by = "cell_id",
      all.x = TRUE
    )
    out[is.na(n_points_in_cell), n_points_in_cell := 0L]
    out[, presence := as.integer(n_points_in_cell > 0L)]
    
    out <- merge(
      out,
      rep_pt,
      by = "cell_id",
      all.x = TRUE
    )
    
    # Ensure lon/lat are NA for absence/background cells
    out[presence == 0L, `:=`(lon = NA_real_, lat = NA_real_)]
    
    sp_dir <- file.path(repo_root, out_root, slug)
    .ensure_dir(sp_dir)
    
    # 1) Full land-grid table
    out_base <- file.path(sp_dir, paste0("occ_", slug, "__grid", cell_km, "km"))
    wrote_full <- FALSE
    if (requireNamespace("arrow", quietly = TRUE)) {
      arrow::write_parquet(as.data.frame(out), paste0(out_base, ".parquet"))
      wrote_full <- TRUE
    } else {
      fwrite(out, paste0(out_base, ".csv"))
      wrote_full <- TRUE
    }
    
    # 2) GEE-ready presence points (ALWAYS write CSV)
    # GEE will create geometries from lon/lat; absences/background handled in EE later.
    gee <- out[presence == 1L & !is.na(lon) & !is.na(lat),
               .(
                 species = sp,
                 lon = as.numeric(lon),
                 lat = as.numeric(lat),
                 cell_id = cell_id,
                 row = as.integer(row),
                 col = as.integer(col),
                 n_points_in_cell = as.integer(n_points_in_cell),
                 cell_km = as.integer(cell_km),
                 crs_grid = crs_grid,
                 policy_tag = policy_tag
               )]
    
    gee_path <- file.path(sp_dir, paste0("presence_points_", cell_km, "km_", slug, ".csv"))
    fwrite(gee, gee_path)
    
    # Optional rasters for quick GIS / sanity checks
    if (terra_ok) {
      r <- terra::rast(
        ncols = ext$ncol, nrows = ext$nrow,
        xmin = ext$xmin, xmax = ext$xmax,
        ymin = ext$ymin, ymax = ext$ymax,
        crs  = crs_grid
      )
      
      ncell <- terra::ncell(r)
      v_presence <- rep(NA_integer_, ncell)
      v_count    <- rep(NA_integer_, ncell)
      
      # terra row index counts from top; our row counts from ymin upward
      terra_row <- ext$nrow - out$row + 1L
      cell_idx  <- terra::cellFromRowCol(r, terra_row, out$col)
      
      v_presence[cell_idx] <- out$presence
      v_count[cell_idx]    <- as.integer(out$n_points_in_cell)
      
      terra::values(r) <- v_presence
      terra::writeRaster(r, file.path(sp_dir, paste0("presence_", cell_km, "km.tif")), overwrite = TRUE)
      
      terra::values(r) <- v_count
      terra::writeRaster(r, file.path(sp_dir, paste0("count_", cell_km, "km.tif")), overwrite = TRUE)
    }
    
    n_points_total <- sum(out$n_points_in_cell)
    n_presence_cells <- sum(out$presence)
    gee_points <- nrow(gee)
    
    if (isTRUE(verbose)) {
      cat("[done] Land cells:", nrow(out),
          " | presence cells:", n_presence_cells,
          " | points total:", n_points_total,
          " | GEE points:", gee_points, "\n")
      cat("[out ] Full grid:", if (grepl("\\.parquet$", list.files(sp_dir, pattern = paste0("^occ_", slug, "__grid", cell_km, "km\\.parquet$"), full.names = FALSE))) "parquet" else "csv", "\n")
      cat("[out ] GEE points:", gee_path, "\n")
    }
    
    summary_dt <- rbind(
      summary_dt,
      data.table(
        species = sp, slug = slug, stage03_file = info$path,
        n_read = n_read,
        n_after_non_na = n_after_non_na,
        n_in_bounds = n_in_bounds,
        n_points_total = as.integer(n_points_total),
        n_presence_cells = as.integer(n_presence_cells),
        gee_points = as.integer(gee_points)
      )
    )
  }
  
  # Summary output
  summary_path_csv <- file.path(repo_root, out_root, "_summary_grid.csv")
  fwrite(summary_dt, summary_path_csv)
  
  if (isTRUE(verbose)) {
    cat("\n============================================================\n")
    cat("Done.\n")
    cat("Outputs:", file.path(repo_root, out_root), "\n")
    cat("Summary:", summary_path_csv, "\n")
    cat("============================================================\n")
  }
  
  invisible(list(
    out_root = file.path(repo_root, out_root),
    grid_land = grid_land,
    summary = summary_dt
  ))
}
