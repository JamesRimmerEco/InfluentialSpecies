# scripts/long_run_wrappers/stage_05_grid_occurrences_all_species.R
#
# ==============================================================================
# Stage 05: Rasterisation / gridding (Stage 04 filtered -> 1 km grid)
# ==============================================================================
#
# Purpose
#   Long-run wrapper to rasterise all species to a 1 km grid (EPSG:3035),
#   using representative observed points (not centroids), with an optional land mask.
#
#   Restart-safe by default:
#     overwrite=FALSE skips species whose outputs already exist.
#
#   IMPORTANT:
#     The Stage 05 engine should backfill summary counts even when skipping
#     (so _summary_grid.csv stays accurate after restart/partial runs).
#
# Engine
#   R/grid_occurrences_stage05_first_observed.R
#
# Inputs
#   data/processed/04_filtered/<slug>/occ_<slug>__filtered.(parquet|rds)
#
# Outputs
#   data/processed/05_grid/<policy_tag>/<slug>/occ_<slug>__grid1km.(parquet|csv)
#   data/processed/05_grid/<policy_tag>/<slug>/presence_points_1km_<slug>.csv
#   data/processed/05_grid/<policy_tag>/_runlog_05_grid.csv
#   data/processed/05_grid/<policy_tag>/_summary_grid.csv
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Find repo root (works from any working directory; also works when pasted) -
find_repo_root <- function(start_dir) {
  markers <- c(".git", "R", "data", "InfluentialSpecies.Rproj", "DESCRIPTION")
  d <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in 1:25) {
    if (any(file.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root <- find_repo_root(getwd())

# ---- Load Stage 05 engine -----------------------------------------------------
engine_path <- file.path(repo_root, "R", "grid_occurrences_stage05_first_observed.R")
if (!file.exists(engine_path)) stop("Stage 05 engine not found: ", engine_path)
source(engine_path)

if (!exists("grid_stage04_to_grid", mode = "function")) {
  stop("Expected function 'grid_stage04_to_grid' not found after sourcing: ", engine_path)
}

# ==============================================================================
# Species list (canonical row set)
# ==============================================================================

species_list_path <- file.path(repo_root, "data", "_meta", "species_list_binomial.csv")
if (!file.exists(species_list_path)) stop("Canonical species list not found: ", species_list_path)

species_names <- readLines(species_list_path, warn = FALSE)
species_names <- species_names[!is.na(species_names)]
species_names <- trimws(gsub("\\s+", " ", species_names))
species_names <- species_names[nzchar(species_names)]
species_names <- species_names[tolower(species_names) != "binomial"]

if (length(species_names) != 100) {
  stop("Canonical species list did not read as 100 names (got ", length(species_names), "). Check: ", species_list_path)
}

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# Inputs/outputs are relative to the repo root
in_root    <- file.path("data", "processed", "04_filtered")
out_stage  <- file.path("data", "processed", "05_grid")

# Output “policy tag” folder name (choose a stable, descriptive label)
policy_tag <- "grid1km_first_observed_landmask_europe_bbox"

# Grid parameters
cell_km   <- 1L
crs_grid  <- "EPSG:3035"

# Europe-ish bbox in lon/lat (WGS84), used as an outer crop window and for context plotting
# (Engine will project internally to the grid CRS.)
bbox_ll <- list(xmin = -12, ymin = 34, xmax = 40, ymax = 72)
# Land mask
use_land_mask <- TRUE

# Optional extras (engine-dependent; safe defaults)
write_geotiff <- FALSE
plot_bbox_map <- TRUE

# Restart behaviour
overwrite <- FALSE

# Logging
write_runlog <- TRUE
verbose <- TRUE

# Optional: limit to first N species for a smoke test
limit_n <- NA_integer_  # e.g. 6L; NA = no limit
if (!is.na(limit_n)) species_names <- species_names[seq_len(min(limit_n, length(species_names)))]

# ==============================================================================
# RUN
# ==============================================================================

grid_stage04_to_grid(
  species_names     = species_names,
  in_root           = in_root,
  out_stage         = out_stage,
  policy_tag        = policy_tag,
  cell_km           = cell_km,
  crs_grid          = crs_grid,
  bbox_ll           = bbox_ll,
  use_land_mask     = use_land_mask,
  write_geotiff     = write_geotiff,
  plot_bbox_map     = plot_bbox_map,
  bbox_map_filename = "bbox_context_map.png",
  overwrite         = overwrite,
  write_runlog      = write_runlog,
  runlog_filename   = "_runlog_05_grid.csv",
  verbose           = verbose
)