# ==============================================================================
# Stage 04 gridding: Stage 03 filtered occurrences -> regular grid (presence/count)
# ==============================================================================
#
# This file contains the "engine" functions for Stage 04.
# A separate wrapper script in /scripts/ sets parameters and calls these functions.
#
# Main entry point:
#   grid_stage03_to_grid()
#
# Optional utilities:
#   plot_europe_bbox_map()
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
      "Tried pairs: ", paste(vapply(coord_pairs, paste, collapse = "/", character(1)), collapse = ", "), "\n",
      "Schema columns include: ", paste(head(schema_names, 60), collapse = ", ")
    )
  }
  
  if (fmt == "rds") {
    dt <- readRDS(path)
    setDT(dt)
    
    if (!("lon" %in% names(dt)) && "decimalLongitude" %in% names(dt)) setnames(dt, "decimalLongitude", "lon")
    if (!("lat" %in% names(dt)) && "decimalLatitude" %in% names(dt)) setnames(dt, "decimalLatitude", "lat")
    if (!("lon" %in% names(dt)) && "longitude" %in% names(dt)) setnames(dt, "longitude", "lon")
    if (!("lat" %in% names(dt)) && "latitude" %in% names(dt)) setnames(dt, "latitude", "lat")
    
    return(dt)
  }
  
  stop("Unknown format: ", fmt)
}

# ---- Project lon/lat -> EPSG:3035 without building sf geometries (fast) -------
.project_ll_to_grid <- function(lon, lat, crs_to) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Need 'sf' for projection: install.packages('sf')")
  }
  
  pts <- cbind(lon, lat)
  
  if ("sf_project" %in% getNamespaceExports("sf")) {
    xy <- sf::sf_project(from = "EPSG:4326", to = crs_to, pts = pts)
    return(list(x = xy[, 1], y = xy[, 2]))
  }
  
  sf_pts <- sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = 4326)
  sf_pts <- sf::st_transform(sf_pts, crs_to)
  xy <- sf::st_coordinates(sf_pts)
  list(x = xy[, 1], y = xy[, 2])
}

# ---- Project EPSG:3035 -> lon/lat (fast when possible) ------------------------
.project_grid_to_ll <- function(x, y, crs_from) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Need 'sf' for projection: install.packages('sf')")
  }
  
  pts <- cbind(x, y)
  
  if ("sf_project" %in% getNamespaceExports("sf")) {
    ll <- sf::sf_project(from = crs_from, to = "EPSG:4326", pts = pts)
    return(list(lon = ll[, 1], lat = ll[, 2]))
  }
  
  sf_pts <- sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = crs_from)
  sf_pts <- sf::st_transform(sf_pts, 4326)
  ll <- sf::st_coordinates(sf_pts)
  list(lon = ll[, 1], lat = ll[, 2])
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

.assign_to_cells <- function(x, y, ext, cell_m) {
  # returns row/col (1-based) and a logical in_bounds
  col <- as.integer(floor((x - ext$xmin) / cell_m) + 1L)
  row <- as.integer(floor((y - ext$ymin) / cell_m) + 1L)
  
  inb <- !is.na(col) & !is.na(row) &
    col >= 1L & col <= ext$ncol &
    row >= 1L & row <= ext$nrow
  
  list(row = row, col = col, in_bounds = inb)
}

# ---- Land polygons (for land masking) ----------------------------------------
.get_land_union <- function(crs_to) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Need 'sf' for land mask: install.packages('sf')")
  }
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
    stop("Need 'rnaturalearth' for land mask: install.packages('rnaturalearth')")
  }
  
  land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  land <- sf::st_make_valid(land)
  land <- sf::st_union(land)
  land <- sf::st_transform(land, crs_to)
  land
}

# ==============================================================================
# Optional: bbox context map
# ==============================================================================

plot_europe_bbox_map <- function(
    bbox_ll,
    out_file,
    pad_deg = 5,
    width_px = 1600,
    height_px = 1200,
    dpi = 200
) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Need 'sf' for bbox map: install.packages('sf')")
  }
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
    stop("Need 'rnaturalearth' for bbox map: install.packages('rnaturalearth')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Need 'ggplot2' for bbox map: install.packages('ggplot2')")
  }
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  world <- sf::st_make_valid(world)
  
  # Bounding box polygon (WGS84 lon/lat)
  bbox_poly <- sf::st_as_sfc(sf::st_bbox(
    c(xmin = bbox_ll$xmin, ymin = bbox_ll$ymin, xmax = bbox_ll$xmax, ymax = bbox_ll$ymax),
    crs = 4326
  ))
  
  # Zoomed-out view limits so the bbox is visible as a rectangle inside the panel
  xlim_view <- c(bbox_ll$xmin - pad_deg, bbox_ll$xmax + pad_deg)
  ylim_view <- c(bbox_ll$ymin - pad_deg, bbox_ll$ymax + pad_deg)
  
  # Breaks (10° major, 5° minor) computed from the view limits
  x_major <- seq(floor(xlim_view[1] / 10) * 10, ceiling(xlim_view[2] / 10) * 10, by = 10)
  x_minor <- seq(floor(xlim_view[1] / 5)  * 5,  ceiling(xlim_view[2] / 5)  * 5,  by = 5)
  y_major <- seq(floor(ylim_view[1] / 10) * 10, ceiling(ylim_view[2] / 10) * 10, by = 10)
  y_minor <- seq(floor(ylim_view[1] / 5)  * 5,  ceiling(ylim_view[2] / 5)  * 5,  by = 5)
  
  # Degree-label formatters
  fmt_lon <- function(x) {
    vapply(x, function(v) {
      if (is.na(v)) return(NA_character_)
      if (v < 0) paste0(abs(v), "\u00B0W") else if (v > 0) paste0(v, "\u00B0E") else "0\u00B0"
    }, character(1))
  }
  
  fmt_lat <- function(y) {
    vapply(y, function(v) {
      if (is.na(v)) return(NA_character_)
      if (v < 0) paste0(abs(v), "\u00B0S") else if (v > 0) paste0(v, "\u00B0N") else "0\u00B0"
    }, character(1))
  }
  
  .ensure_dir(dirname(out_file))
  
  gg <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = world, fill = NA, linewidth = 0.3) +
    # BBox drawn as dashed so it reads as an overlay
    ggplot2::geom_sf(data = bbox_poly, fill = NA, linewidth = 1.0, linetype = "dashed") +
    # Turn off coord_sf graticule (the source of the two dark lines)
    # ...and use explicit scales below for ticks + gridlines.
    ggplot2::coord_sf(
      xlim = xlim_view,
      ylim = ylim_view,
      expand = FALSE,
      datum = NA
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_major,
      minor_breaks = x_minor,
      labels = fmt_lon
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_major,
      minor_breaks = y_minor,
      labels = fmt_lat
    ) +
    ggplot2::labs(
      title = "Europe context map with modelling bounding box",
      x = "Longitude", y = "Latitude"
    )
  
  ggplot2::ggsave(
    filename = out_file,
    plot = gg,
    width = width_px / dpi,
    height = height_px / dpi,
    dpi = dpi
  )
  
  invisible(out_file)
}

# ==============================================================================
# Main engine: Stage 03 -> grid
# ==============================================================================

grid_stage03_to_grid <- function(
    species_names,
    in_root = file.path("data", "processed", "03_filtered"),
    out_stage = "04_grid",
    policy_tag = "grid25km_test_landmask_europe_bbox",
    cell_km = 25,
    crs_grid = "EPSG:3035",
    bbox_ll = list(xmin = -25, xmax = 45, ymin = 34, ymax = 72),
    use_land_mask = TRUE,
    write_geotiff = TRUE,
    plot_bbox_map = FALSE,
    bbox_map_filename = "bbox_context_map.png",
    verbose = TRUE
) {
  repo_root <- normalizePath(getwd())
  
  cell_m <- cell_km * 1000
  out_root <- file.path("data", "processed", out_stage, policy_tag)
  .ensure_dir(file.path(repo_root, out_root))
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("This stage requires 'sf'. Install it first: install.packages('sf')")
  }
  
  # ---------------------------------------------------------------------------
  # Optional bbox context map (lon/lat)
  # ---------------------------------------------------------------------------
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
    cat("============================================================\n")
  }
  
  # Save the land grid once (useful for debugging and later reuse)
  grid_out_base <- file.path(repo_root, out_root, paste0("grid", cell_km, "km_land_cells"))
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(as.data.frame(grid_land), paste0(grid_out_base, ".parquet"))
  } else {
    fwrite(grid_land, paste0(grid_out_base, ".csv"))
  }
  
  # ---------------------------------------------------------------------------
  # Per-species gridding
  # ---------------------------------------------------------------------------
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
    n_presence_cells = integer()
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
          n_points_total = 0L, n_presence_cells = 0L
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
          n_points_total = 0L, n_presence_cells = 0L
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
          n_points_total = 0L, n_presence_cells = 0L
        )
      )
      next
    }
    
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
          n_points_total = 0L, n_presence_cells = 0L
        )
      )
      next
    }
    
    counts <- dt[, .(n_points_in_cell = .N), by = .(row, col)]
    counts[, cell_id := sprintf("r%05d_c%05d", row, col)]
    
    out <- merge(
      grid_land[, .(cell_id, row, col, x_center, y_center, on_land)],
      counts[, .(cell_id, n_points_in_cell)],
      by = "cell_id",
      all.x = TRUE
    )
    out[is.na(n_points_in_cell), n_points_in_cell := 0L]
    out[, presence := as.integer(n_points_in_cell > 0L)]
    
    ll <- .project_grid_to_ll(out$x_center, out$y_center, crs_from = crs_grid)
    out[, `:=`(
      lon_center = ll$lon,
      lat_center = ll$lat
    )]
    
    sp_dir <- file.path(repo_root, out_root, slug)
    .ensure_dir(sp_dir)
    
    out_base <- file.path(sp_dir, paste0("occ_", slug, "__grid", cell_km, "km"))
    if (requireNamespace("arrow", quietly = TRUE)) {
      arrow::write_parquet(as.data.frame(out), paste0(out_base, ".parquet"))
    } else {
      fwrite(out, paste0(out_base, ".csv"))
    }
    
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
    
    if (isTRUE(verbose)) {
      cat("[done] Land cells:", nrow(out),
          " | presence cells:", n_presence_cells,
          " | points total:", n_points_total, "\n")
    }
    
    summary_dt <- rbind(
      summary_dt,
      data.table(
        species = sp, slug = slug, stage03_file = info$path,
        n_read = n_read,
        n_after_non_na = n_after_non_na,
        n_in_bounds = n_in_bounds,
        n_points_total = as.integer(n_points_total),
        n_presence_cells = as.integer(n_presence_cells)
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Summary output
  # ---------------------------------------------------------------------------
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
