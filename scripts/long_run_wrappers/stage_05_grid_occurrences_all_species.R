# scripts/long_run_wrappers/stage_05_grid_occurrences_all_species.R
#
# ==============================================================================
# Stage 05: Rasterisation / gridding (Stage 04 filtered -> 1 km grid)
# ==============================================================================
#
# Purpose:
#   Long-run “home run” wrapper to rasterise all species to a 1 km grid (EPSG:3035),
#   using representative observed points (not centroids), with an optional land mask.
#
#   Restart-safe by default: overwrite=FALSE will skip species with existing outputs.
#
# Engine:
#   R/grid_occurrences_stage04_first_observed.R
#
# Inputs:
#   data/processed/04_filtered/<slug>/occ_<slug>__filtered.(parquet|rds)
#
# Outputs:
#   data/processed/05_grid/<policy_tag>/<slug>/occ_<slug>__grid1km.(parquet|csv)
#   data/processed/05_grid/<policy_tag>/<slug>/presence_points_1km_<slug>.csv
#   data/processed/05_grid/<policy_tag>/_runlog_05_grid.csv
#   data/processed/05_grid/<policy_tag>/_summary_grid.csv
#
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Find repo root (works from any scripts/ subfolder) -----------------------
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file) || !nzchar(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run via source('.../scripts/.../stage_05_grid_occurrences_all_species.R') from a file, not copy/paste."
  )
}

script_dir <- dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE))

find_repo_root <- function(start_dir) {
  markers <- c(".git", "R", "data", "InfluentialSpecies.Rproj", "DESCRIPTION")
  d <- start_dir
  for (i in 1:20) {
    if (any(file.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root <- find_repo_root(script_dir)

# ---- Load the Stage 05 gridding engine (deterministic) ------------------------
engine_fn <- file.path(repo_root, "R", "grid_occurrences_stage05_first_observed.R")
if (!file.exists(engine_fn)) {
  stop("Can't find Stage 05 gridding engine at: ", engine_fn)
}
source(engine_fn)

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# ---- Species list -------------------------------------------------------------
species_csv <- file.path(repo_root, "data", "_meta", "species_list_binomial.csv")
if (!file.exists(species_csv)) stop("Can't find species list at: ", species_csv)

species_names <- read.csv(species_csv, stringsAsFactors = FALSE, header = TRUE)[[1]]
species_names <- as.character(species_names)
species_names <- trimws(species_names)
species_names <- species_names[!is.na(species_names) & nzchar(species_names)]
species_names <- species_names[tolower(species_names) != "binomial"]  # belt + braces
species_names <- unique(species_names)

if (length(species_names) < 10) {
  stop("Species list looks unexpectedly short (", length(species_names), "). Check: ", species_csv)
}

# ---- Stage handoff ------------------------------------------------------------
# Stage 04 -> Stage 05
in_root <- file.path("data", "processed", "04_filtered")

# ---- Output identity ----------------------------------------------------------
cell_km <- 1
use_land_mask <- TRUE
out_stage <- "05_grid"

policy_tag <- paste0("grid", cell_km, "km_first_observed_",
                     if (use_land_mask) "landmask_" else "",
                     "europe_bbox")

# ---- Grid settings ------------------------------------------------------------
crs_grid <- "EPSG:3035"  # ETRS89 / LAEA Europe

bbox_ll <- list(
  xmin = -25,  # lon
  xmax =  45,
  ymin =  34,  # lat
  ymax =  72
)

# Optional outputs
write_geotiff <- FALSE        # 1 km GeoTIFFs can be large; enable only if needed
plot_bbox_map <- TRUE         # writes a single context map per run (policy_tag folder)

# Long-run behaviour
overwrite <- FALSE            # restart-safe default

write_runlog <- TRUE          # writes _runlog_05_grid.csv incrementally
verbose <- TRUE

# Optional: limit to first N species for a smoke test
limit_n <- NA_integer_        # e.g. 6L; NA = no limit
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