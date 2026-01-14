# ==============================================================================
# Script: 06_prepare_gee_export.R
# Purpose: Prepare spatial data for Google Earth Engine integration
# Author: [Author name]
# Date: [Date]
# ==============================================================================

# Description:
# This script formats spatial occurrence data and analysis outputs for upload
# to Google Earth Engine (GEE). It handles coordinate precision, file format
# conversion, size optimization, and metadata documentation required for
# GEE asset management.

# ==============================================================================
# Setup
# ==============================================================================

# Load required libraries
library(tidyverse)
library(sf)
library(here)
library(glue)
library(jsonlite)

# Source project configuration
source(here("scripts", "02_configure_project.R"))

# ==============================================================================
# Load Data
# ==============================================================================

message("Loading spatial data for GEE export...")

# Load processed spatial data
# spatial_data <- st_read(here(paths$processed, "occurrences_spatial.geojson"))
# grid_data <- st_read(here(paths$processed, "spatial_grid_aggregated.geojson"))

# ==============================================================================
# Data Formatting for GEE
# ==============================================================================

#' Format coordinates for GEE precision requirements
#'
#' @param sf_obj sf object
#' @param precision Number of decimal places for coordinates
#' @return sf object with rounded coordinates
#'
format_coordinates <- function(sf_obj, precision = 6) {
  
  message(glue("Formatting coordinates to {precision} decimal places..."))
  
  # Extract coordinates, round, and reconstruct geometry
  coords <- st_coordinates(sf_obj)
  coords_rounded <- round(coords, precision)
  
  # Create new geometry with rounded coordinates
  # Note: This is a simplified approach; actual implementation depends on geometry type
  sf_rounded <- sf_obj
  st_geometry(sf_rounded) <- st_sfc(
    st_multipoint(coords_rounded),
    crs = st_crs(sf_obj)
  )
  
  return(sf_rounded)
}

#' Simplify attribute table for GEE
#'
#' @param sf_obj sf object
#' @return sf object with simplified attributes
#'
simplify_attributes <- function(sf_obj) {
  
  message("Simplifying attribute table...")
  
  # Keep only essential columns
  # GEE has attribute limitations, so reduce columns to necessary fields
  # essential_cols <- c("species", "date", "record_id", "geometry")
  # 
  # sf_simplified <- sf_obj %>%
  #   select(any_of(essential_cols))
  
  return(sf_obj)
}

#' Convert dates to GEE-compatible format
#'
#' @param sf_obj sf object with date column
#' @return sf object with formatted dates
#'
format_dates_for_gee <- function(sf_obj) {
  
  message("Converting dates to GEE-compatible format...")
  
  # GEE uses milliseconds since Unix epoch
  # sf_formatted <- sf_obj %>%
  #   mutate(
  #     date_gee = as.numeric(as.POSIXct(date)) * 1000
  #   )
  
  return(sf_obj)
}

# ==============================================================================
# GEE Asset Preparation
# ==============================================================================

#' Prepare point features for GEE FeatureCollection
#'
#' @param sf_obj sf point object
#' @param asset_name Name for the GEE asset
#' @return Formatted sf object and metadata
#'
prepare_point_asset <- function(sf_obj, asset_name) {
  
  message(glue("Preparing point asset: {asset_name}"))
  
  # Format data
  formatted <- sf_obj %>%
    format_coordinates(precision = options$gee_precision) %>%
    simplify_attributes() %>%
    format_dates_for_gee()
  
  # Create metadata
  metadata <- list(
    asset_name = asset_name,
    asset_type = "FeatureCollection",
    geometry_type = "Point",
    n_features = nrow(formatted),
    crs = st_crs(formatted)$input,
    date_created = Sys.time(),
    source_script = "06_prepare_gee_export.R",
    properties = names(formatted)
  )
  
  return(list(data = formatted, metadata = metadata))
}

#' Prepare polygon features for GEE FeatureCollection
#'
#' @param sf_obj sf polygon object
#' @param asset_name Name for the GEE asset
#' @return Formatted sf object and metadata
#'
prepare_polygon_asset <- function(sf_obj, asset_name) {
  
  message(glue("Preparing polygon asset: {asset_name}"))
  
  # Format data
  formatted <- sf_obj %>%
    format_coordinates(precision = options$gee_precision) %>%
    simplify_attributes()
  
  # Create metadata
  metadata <- list(
    asset_name = asset_name,
    asset_type = "FeatureCollection",
    geometry_type = "Polygon",
    n_features = nrow(formatted),
    crs = st_crs(formatted)$input,
    date_created = Sys.time(),
    source_script = "06_prepare_gee_export.R",
    properties = names(formatted)
  )
  
  return(list(data = formatted, metadata = metadata))
}

# ==============================================================================
# Size Optimization
# ==============================================================================

#' Check file size and warn if too large for GEE
#'
#' @param file_path Path to file
#' @param max_size_mb Maximum size in MB
#'
check_file_size <- function(file_path, max_size_mb = 250) {
  
  size_mb <- file.size(file_path) / (1024^2)
  
  message(glue("File size: {round(size_mb, 2)} MB"))
  
  if (size_mb > max_size_mb) {
    warning(glue("File size ({round(size_mb, 2)} MB) exceeds recommended maximum ({max_size_mb} MB)."))
    message("Consider:")
    message("  - Splitting into multiple assets")
    message("  - Reducing attribute columns")
    message("  - Aggregating to coarser spatial resolution")
  }
  
  return(size_mb)
}

# ==============================================================================
# Export Functions
# ==============================================================================

#' Export asset to GeoJSON format
#'
#' @param asset_list List with data and metadata
#' @param output_dir Output directory
#'
export_geojson_asset <- function(asset_list, output_dir) {
  
  asset_name <- asset_list$metadata$asset_name
  
  # Save GeoJSON
  geojson_path <- here(output_dir, glue("{asset_name}.geojson"))
  st_write(
    asset_list$data, 
    geojson_path, 
    delete_dsn = TRUE,
    quiet = TRUE
  )
  
  message(glue("Exported: {geojson_path}"))
  
  # Check file size
  size <- check_file_size(geojson_path)
  
  # Save metadata
  metadata_path <- here(output_dir, glue("{asset_name}_metadata.json"))
  write_json(asset_list$metadata, metadata_path, pretty = TRUE, auto_unbox = TRUE)
  
  message(glue("Metadata: {metadata_path}"))
  
  return(geojson_path)
}

#' Export asset to Shapefile format (alternative)
#'
#' @param asset_list List with data and metadata
#' @param output_dir Output directory
#'
export_shapefile_asset <- function(asset_list, output_dir) {
  
  asset_name <- asset_list$metadata$asset_name
  
  # Save Shapefile
  shp_path <- here(output_dir, glue("{asset_name}.shp"))
  st_write(
    asset_list$data,
    shp_path,
    delete_dsn = TRUE,
    quiet = TRUE
  )
  
  message(glue("Exported: {shp_path}"))
  
  # Save metadata
  metadata_path <- here(output_dir, glue("{asset_name}_metadata.json"))
  write_json(asset_list$metadata, metadata_path, pretty = TRUE, auto_unbox = TRUE)
  
  return(shp_path)
}

# ==============================================================================
# Execute GEE Export Preparation
# ==============================================================================

message("\n=== Preparing Data for Google Earth Engine ===\n")

# Prepare occurrence points asset
# occurrences_asset <- prepare_point_asset(
#   spatial_data,
#   asset_name = "influential_species_occurrences"
# )

# Prepare spatial grid asset
# grid_asset <- prepare_polygon_asset(
#   grid_data,
#   asset_name = "influential_species_grid"
# )

# ==============================================================================
# Export Assets
# ==============================================================================

# Export to GeoJSON (recommended for GEE)
# if (options$gee_format == "GeoJSON") {
#   export_geojson_asset(occurrences_asset, paths$gee_exports)
#   export_geojson_asset(grid_asset, paths$gee_exports)
# }

# Export to Shapefile (alternative)
# if (options$gee_format == "Shapefile") {
#   export_shapefile_asset(occurrences_asset, paths$gee_exports)
#   export_shapefile_asset(grid_asset, paths$gee_exports)
# }

# ==============================================================================
# Generate Upload Instructions
# ==============================================================================

#' Create GEE upload instructions document
#'
create_upload_instructions <- function() {
  
  instructions <- glue("
# Google Earth Engine Upload Instructions

## Generated: {Sys.time()}

## Files Ready for Upload
Located in: {paths$gee_exports}

## Upload Process

### Option 1: GEE Code Editor Asset Manager
1. Open Google Earth Engine Code Editor (https://code.earthengine.google.com/)
2. Navigate to the Assets tab
3. Click 'NEW' > 'Table Upload' (for vector data)
4. Select the .geojson or .shp file
5. Configure asset properties:
   - Asset ID: Choose a meaningful name
   - Properties: Review and confirm attribute columns
6. Click 'Upload'

### Option 2: GEE Python API
```python
import ee
ee.Initialize()

# Upload asset
task = ee.batch.Export.table.toAsset(
  collection=ee.FeatureCollection('path/to/local/file.geojson'),
  description='Upload influential species data',
  assetId='users/YOUR_USERNAME/influential_species_occurrences'
)
task.start()
```

### Option 3: Command Line (earthengine CLI)
```bash
earthengine upload table --asset_id=users/YOUR_USERNAME/asset_name path/to/file.geojson
```

## Asset Naming Conventions
- Use descriptive, lowercase names with underscores
- Include date or version if applicable
- Example: influential_species_occurrences_2024

## After Upload
1. Check asset properties in Code Editor
2. Verify geometry and attributes loaded correctly
3. Update GEE script to reference new asset ID
4. Document asset ID in project documentation

## Troubleshooting
- File size limit: 250 MB per asset
- If file too large, consider splitting by species or region
- Check coordinate system is WGS84 (EPSG:4326)
- Ensure attribute names are GEE-compatible (no special characters)

For more information: https://developers.google.com/earth-engine/guides/asset_manager
")
  
  instructions_file <- here(paths$gee_exports, "GEE_UPLOAD_INSTRUCTIONS.md")
  writeLines(instructions, instructions_file)
  
  message(glue("Upload instructions saved to: {instructions_file}"))
}

# Generate instructions
create_upload_instructions()

message("\n=== Script Notes ===")
message("This is an annotated placeholder script.")
message("Review exported files before uploading to GEE.")
message("Ensure compliance with GEE terms of service and data policies.")
message("Document GEE asset IDs for use in subsequent workflows.")
message("See docs/gee_integration.md for detailed guidance.")

# ==============================================================================
# End of Script
# ==============================================================================
