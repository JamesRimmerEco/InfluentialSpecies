# ==============================================================================
# Script: 02_configure_project.R
# Purpose: Set project parameters, paths, and configuration options
# Author: [Author name]
# Date: [Date]
# ==============================================================================

# Description:
# This script defines project-wide parameters, file paths, and configuration
# settings. Modify these parameters to suit your specific research question
# and data sources. This script should be sourced at the beginning of all
# analysis scripts.

# ==============================================================================
# Load Required Libraries
# ==============================================================================

library(here)
library(glue)

# ==============================================================================
# Project Paths
# ==============================================================================

# Define directory structure using 'here' for reproducible paths
paths <- list(
  # Data directories
  raw_data      = here("data", "raw"),
  processed     = here("data", "processed"),
  outputs       = here("data", "outputs"),
  
  # Script and function directories
  scripts       = here("scripts"),
  functions     = here("R"),
  
  # Documentation
  docs          = here("docs"),
  
  # Specific output subdirectories
  figures       = here("data", "outputs", "figures"),
  tables        = here("data", "outputs", "tables"),
  gee_exports   = here("data", "outputs", "gee_exports")
)

# Create directories if they don't exist
for (path in paths) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created directory: ", path)
  }
}

# ==============================================================================
# Project Parameters
# ==============================================================================

# Study region parameters
# NOTE: These are placeholders - adjust based on your research question
study_parameters <- list(
  # Geographic extent (example: UK coordinates in WGS84)
  # Modify these bounds to match your study area
  bbox = list(
    xmin = -10.0,  # Western longitude
    xmax = 2.0,    # Eastern longitude
    ymin = 49.0,   # Southern latitude
    ymax = 61.0    # Northern latitude
  ),
  
  # Coordinate reference system
  crs_wgs84 = "EPSG:4326",      # WGS84 for GBIF data
  crs_projected = "EPSG:27700",  # British National Grid (example)
  
  # Temporal extent
  # Adjust dates based on your research requirements
  date_start = "2000-01-01",
  date_end = "2023-12-31"
)

# Species of interest
# NOTE: This is a placeholder list - replace with your target species
# Use scientific names that match the taxonomic backbone you'll query
species_list <- c(
  # Example species - replace with your study organisms
  "Species name 1",
  "Species name 2",
  "Species name 3"
)

# Data quality parameters
quality_parameters <- list(
  # Coordinate uncertainty threshold (metres)
  max_coord_uncertainty = 1000,
  
  # Minimum number of records per species
  min_records = 10,
  
  # Basis of record to include (GBIF terms)
  acceptable_basis = c(
    "HUMAN_OBSERVATION",
    "OBSERVATION",
    "MACHINE_OBSERVATION",
    "PRESERVED_SPECIMEN"
  )
)

# ==============================================================================
# Analysis Options
# ==============================================================================

# Set random seed for reproducibility
set.seed(42)

# Processing options
options <- list(
  # Parallel processing
  use_parallel = FALSE,  # Set to TRUE if using parallel computation
  n_cores = 4,           # Number of cores for parallel processing
  
  # Output formats
  save_rds = TRUE,       # Save R data objects
  save_csv = TRUE,       # Save CSV files
  save_geojson = TRUE,   # Save spatial data as GeoJSON
  
  # GEE export settings
  gee_format = "GeoJSON", # Format for GEE assets (GeoJSON, Shapefile)
  gee_precision = 6       # Decimal places for coordinates
)

# ==============================================================================
# File Naming Conventions
# ==============================================================================

# Generate standardised file names with timestamps
create_filename <- function(prefix, extension, add_timestamp = TRUE) {
  if (add_timestamp) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- glue("{prefix}_{timestamp}.{extension}")
  } else {
    filename <- glue("{prefix}.{extension}")
  }
  return(filename)
}

# ==============================================================================
# Export Configuration
# ==============================================================================

# Save configuration to file for documentation
config_export <- list(
  study_parameters = study_parameters,
  species_list = species_list,
  quality_parameters = quality_parameters,
  options = options,
  paths = paths,
  timestamp = Sys.time()
)

saveRDS(
  config_export,
  file = here("data", "processed", "project_config.rds")
)

message("Project configuration complete!")
message("Configuration saved to: ", here("data", "processed", "project_config.rds"))

# ==============================================================================
# End of Script
# ==============================================================================
