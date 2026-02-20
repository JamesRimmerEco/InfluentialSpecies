# ==============================================================================
# FILE: scripts/test_wrappers/grid_stage03_25km_test_v02_first_observed.R
# ==============================================================================
#
# ==============================================================================
# Rasterisation / gridding: Stage 03 filtered occurrences -> regular grid
# ==============================================================================
#
# Wrapper for the Stage 04 gridding engine.
#
# What you change here:
#   - which species to grid (test subset vs full list)
#   - cell_km resolution (25 km now; later experiments can change this)
#   - bbox extent and land mask toggle
#   - output policy_tag (keeps outputs organised and prevents accidental mixing)
#
# This wrapper calls the engine:
#   R/grid_occurrences_stage03_v02_first_observed.R
#
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- Find repo root (works from any scripts/ subfolder) ----
this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(this_file) || !nzchar(this_file)) {
  stop(
    "Can't determine script path (sys.frame(1)$ofile is NULL). ",
    "Run via source('.../scripts/.../grid_stage03_25km_test_v02_first_observed.R') from a file, not copy/paste."
  )
}
script_dir <- dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE))

find_repo_root <- function(start_dir) {
  markers <- c(".git", "R", "data", "InfluentialSpecies.Rproj", "DESCRIPTION")
  d <- start_dir
  for (i in 1:15) {
    if (any(file.exists(file.path(d, markers)))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Couldn't find repo root walking up from: ", start_dir)
}

repo_root <- find_repo_root(script_dir)

# ==============================================================================
# CONTROL PANEL
# ==============================================================================

# ---- Species list -------------------------------------------------------------
# For quick tests you can hard-code a small subset.
# For a full run, read the canonical binomial list produced in meta.
use_meta_species_list <- FALSE # Change to TRUE for real pull

meta_species_csv <- file.path(repo_root, "data", "_meta", "species_list_binomial.csv")

if (isTRUE(use_meta_species_list)) {
  if (!requireNamespace("readr", quietly = TRUE)) stop("This wrapper needs readr: install.packages('readr')")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("This wrapper needs stringr: install.packages('stringr')")
  
  x <- readr::read_csv(meta_species_csv, show_col_types = FALSE)[[1]]
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- x[!is.na(x) & nzchar(x)]
  x <- x[tolower(x) != "binomial"]  # belt + braces
  x <- x[!duplicated(x)]
  
  # Fail fast if the meta file is wrong (keeps long runs from silently doing nothing)
  ok_binom <- grepl("^[A-Z][a-z-]+\\s+[a-z-]+$", x)
  if (!all(ok_binom)) {
    bad <- x[!ok_binom]
    stop(
      "Meta species list contains non-binomials. Examples: ",
      paste(utils::head(bad, 10), collapse = " | "),
      "\nCheck: ", meta_species_csv
    )
  }
  
  species_names <- x
} else {
  species_names <- c("Andrena fulva")
}


# Optional: limit to first N species for a smoke test
limit_n <- NA_integer_  # e.g. 6L; NA = no limit
if (!is.na(limit_n)) species_names <- species_names[seq_len(min(limit_n, length(species_names)))]

# ---- Where to read Stage 03 from ---------------------------------------------
in_root <- file.path("data", "processed", "03_filtered")

# ---- Where to write gridded outputs ------------------------------------------
out_stage <- "04_grid"

# Policy tag should capture the important choices that define output identity.
# This version uses representative observed points (NOT centroids).
cell_km  <- 25
use_land_mask <- TRUE
policy_tag <- paste0("grid", cell_km, "km_first_observed_landmask_europe_bbox")

# ---- Grid settings ------------------------------------------------------------
crs_grid <- "EPSG:3035"  # ETRS89 / LAEA Europe

# Europe-ish extent as a bounding box in lon/lat (WGS84).
bbox_ll <- list(
  xmin = -25,  # lon
  xmax =  45,
  ymin =  34,  # lat
  ymax =  72
)

# Optional outputs
write_geotiff <- TRUE
plot_bbox_map <- TRUE

# ==============================================================================
# RUN
# ==============================================================================

source(file.path(repo_root, "R", "grid_occurrences_stage03_v02_first_observed.R"))

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
