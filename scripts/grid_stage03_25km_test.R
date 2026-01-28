# ==============================================================================
# Rasterisation / gridding: what this script is doing 
# ==============================================================================
#
# Goal
# ----
# We currently have a big table of occurrence points for each species (one row = one observation
# with a longitude/latitude). For presence/absence modelling we don’t want millions of points.
# We want a consistent grid over space where each grid cell is either:
#   - presence = 1  (the species was recorded in that cell at least once)
#   - presence = 0  (the species was not recorded in that cell)
#
# Important: "absence" here does NOT mean the species is truly absent.
# It just means we do not have a record in that grid cell.
#
#
# Why make a grid at all?
# -----------------------
# 1) It reduces duplication and sampling bias. Thousands of points in the same town become
#    one “presence cell” instead of overwhelming the model.
# 2) It gives us a consistent unit (e.g., 25 km x 25 km) across Europe, so species can be compared.
# 3) It creates an easy format for later steps (MaxEnt input, GIS, AlphaEarth-style pipelines).
#
#
# Why 25 km x 25 km? 
# ------------------
# Because our study area is Europe-wide and we have lots of records for some species.
# 25 km cells are a  “first pass” resolution, but this is arbitrary and might change. 
#
#
# Why we switch coordinate systems (CRS - coordinate reference system)
# --------------------------------------
# Longitude/latitude degrees are NOT a fixed distance (a degree of longitude gets smaller
# as you go north). If we tried to make “25 km cells” in degrees we would get uneven cells.
#
# So we project the data into a Europe-wide metric coordinate system:
#   EPSG:3035 (ETRS89 / LAEA Europe)
# In this CRS, distances are in metres, so “25 km” is literally 25,000 metres everywhere.
#
#
# What the grid actually is
# -------------------------
# The grid is just a tiled set of squares covering the study area. Each square has:
#   - a row number and column number (grid indices)
#   - a centre coordinate (x_center, y_center) in metres (EPSG:3035)
#   - an ID like "r00012_c00034"
#
# We build the grid once, then reuse it for every species.
#
#
# What “Europe-ish bounding box” means here
# -----------------------------------------
# This script needs a study area to build the grid over.
# For now we use a simple bounding box in lon/lat (WGS84), roughly covering Europe.
#
# This bounding box is not a scientific definition of Europe.
# It’s simply a practical extent so the pipeline runs end-to-end.
# We can tighten/adjust it later if needed.
#
#
# Why we use a land mask
# --------------------------------------------------
# Many grid cells in the Europe bounding box will be ocean.
# For a terrestrial species model, those cells are unhelpful clutter.
#
# The “land mask” step removes ocean cells using a simple and fast rule:
#   keep a cell if the CENTRE POINT of the cell falls on land polygons.
#
# This is a heuristic:
#   - it will correctly remove obvious ocean
#   - it might keep/remove some coastal edge cells imperfectly
# For a first pass, that’s fine.
#
#
# How we turn points into grid cells
# ---------------------------------
# For each species:
#
# 1) Read the Stage 03 filtered dataset for that species.
#    This is already cleaned and policy-filtered.
#
# 2) Keep only coordinates (lon/lat) and drop any missing coordinates.
#
# 3) Project lon/lat into EPSG:3035 (so they are in metres).
#
# 4) For each point, work out which grid cell it falls in.
#    We do this by calculating the row/col index from the x/y position.
#
# 5) Count how many points fall into each cell:
#      n_points_in_cell
#
# 6) Convert counts into presence/absence:
#      presence = 1 if n_points_in_cell > 0
#      presence = 0 if n_points_in_cell == 0
#
# 7) Join those results onto the full LAND grid so we get:
#    - all land cells (background space) plus
#    - presence and count information
#
#
# What we output
# ------------------------
# The core output is a grid table for each species (one row per land cell):
#   cell_id, row, col, x_center, y_center, lon_center, lat_center,
#   n_points_in_cell, presence
#
# Optionally, if the terra package is installed, we also write GeoTIFF rasters:
#   - presence_25km.tif  (1/0/NA)
#   - count_25km.tif     (counts/NA)
#
# In the rasters:
#   - NA = not in our modelling grid (e.g. sea cells, if land mask is on)
#   - 0  = land cell, but no record for this species
#   - 1  = land cell with at least one record
#
# What “presence = 0” does and does not mean
# ------------------------------------------
# presence = 0 means “no observation in this cell in our dataset”.
# It does NOT mean the species is truly absent there.
# That distinction matters: this is still presence-only data,
# we are just converting it into a grid representation.
#
#
# Quick sanity checks after running
# --------------------------------------------
# For each species:
#   - How many land cells are there total?
#   - How many cells have presence = 1?
#   - Does sum(n_points_in_cell) roughly match the number of Stage 03 points?
#   - Plot a quick map (presence cells should look geographically sensible)
#
# If presence cells look wrong (e.g. lots of points in the sea), then either:
#   - coordinates are wrong (upstream issue), or
#   - projection / grid extent is wrong (this script), or
#   - land mask is too permissive/strict (coastline edge effects)
#
# For now, the goal is just a clean, reproducible grid output for testing.
# We can refine “what is the correct study extent / background area?” later.
# ==============================================================================


suppressPackageStartupMessages({
  library(data.table)
})

# ---- Find this script’s directory (works when sourced from a file) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run this via source('.../scripts/grid_stage03_25km_test.R') from a file, not by copy/paste."
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
# Stage folder is generic ("04_grid"). Resolution/choices sit in the policy_tag subfolder.
out_stage <- "04_grid"
policy_tag <- "grid25km_test_landmask_europe_bbox"

# ---- Grid settings ------------------------------------------------------------
cell_km  <- 25
crs_grid <- "EPSG:3035"  # ETRS89 / LAEA Europe

# Europe-ish extent as a bounding box in lon/lat (WGS84).
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

# Optional: write a context map of Europe with the bbox drawn (for speed make FALSE)
plot_bbox_map <- TRUE

# ==============================================================================
# RUN
# ==============================================================================

source(file.path("R", "grid_occurrences_stage03.R"))

grid_stage03_to_grid(
  species_names   = species_names,
  in_root         = in_root,
  out_stage       = out_stage,
  policy_tag      = policy_tag,
  cell_km         = cell_km,
  crs_grid        = crs_grid,
  bbox_ll         = bbox_ll,
  use_land_mask   = use_land_mask,
  write_geotiff   = write_geotiff,
  plot_bbox_map   = plot_bbox_map,
  bbox_map_filename = "bbox_context_map.png",
  verbose         = TRUE
)
