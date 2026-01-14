# ==============================================================================
# Script: 05_spatial_processing.R
# Purpose: Conduct spatial analyses and generate spatial summaries
# Author: [Author name]
# Date: [Date]
# ==============================================================================

# Description:
# This script performs spatial processing on cleaned occurrence data, including
# coordinate transformations, spatial aggregation, grid-based summaries, and
# calculation of spatial metrics. Outputs are designed to be compatible with
# further analysis or GEE integration.

# ==============================================================================
# Setup
# ==============================================================================

# Load required libraries
library(tidyverse)
library(sf)
library(terra)
library(here)
library(glue)
library(tmap)

# Source project configuration
source(here("scripts", "02_configure_project.R"))

# ==============================================================================
# Load Data
# ==============================================================================

message("Loading cleaned occurrence data...")

# Load cleaned data from previous script
# clean_data <- readRDS(here(paths$processed, "occurrences_clean_[timestamp].rds"))

# ==============================================================================
# Convert to Spatial Objects
# ==============================================================================

#' Convert occurrence dataframe to sf object
#'
#' @param data Dataframe with longitude and latitude columns
#' @param crs Coordinate reference system (EPSG code)
#' @return sf object
#'
convert_to_sf <- function(data, crs = study_parameters$crs_wgs84) {
  
  message("Converting to spatial features...")
  
  # Create sf object
  data_sf <- st_as_sf(
    data,
    coords = c("longitude", "latitude"),
    crs = crs,
    remove = FALSE
  )
  
  message(glue("Created sf object with {nrow(data_sf)} features"))
  
  return(data_sf)
}

#' Project spatial data to different CRS
#'
project_spatial_data <- function(sf_obj, target_crs) {
  
  message(glue("Projecting to CRS: {target_crs}"))
  
  projected <- st_transform(sf_obj, crs = target_crs)
  
  return(projected)
}

# ==============================================================================
# Spatial Aggregation
# ==============================================================================

#' Create spatial grid for aggregation
#'
#' @param bbox Bounding box (xmin, ymin, xmax, ymax)
#' @param cell_size Grid cell size in map units
#' @param crs Coordinate reference system
#' @return sf object of grid polygons
#'
create_spatial_grid <- function(bbox, cell_size, crs) {
  
  message(glue("Creating spatial grid with cell size: {cell_size}"))
  
  # NOTE: Implement grid creation
  # Example approach using sf:
  # grid <- st_make_grid(
  #   x = st_bbox(c(xmin = bbox$xmin, ymin = bbox$ymin, 
  #                 xmax = bbox$xmax, ymax = bbox$ymax), crs = crs),
  #   cellsize = cell_size,
  #   what = "polygons"
  # ) %>%
  #   st_as_sf() %>%
  #   mutate(grid_id = row_number())
  
  # return(grid)
  
  return(NULL)
}

#' Aggregate occurrences to grid cells
#'
#' @param occurrences sf object of occurrence points
#' @param grid sf object of grid polygons
#' @return sf object with aggregated counts per grid cell
#'
aggregate_to_grid <- function(occurrences, grid) {
  
  message("Aggregating occurrences to grid cells...")
  
  # Spatial join occurrences to grid
  # aggregated <- st_join(grid, occurrences) %>%
  #   group_by(grid_id) %>%
  #   summarise(
  #     n_records = n(),
  #     n_species = n_distinct(species),
  #     species_list = paste(unique(species), collapse = "; "),
  #     .groups = "drop"
  #   )
  
  # return(aggregated)
  
  return(NULL)
}

# ==============================================================================
# Spatial Metrics
# ==============================================================================

#' Calculate species range metrics
#'
#' @param occurrences sf object of occurrence points
#' @param species_name Target species
#' @return List of spatial metrics
#'
calculate_range_metrics <- function(occurrences, species_name) {
  
  # Filter to target species
  species_data <- occurrences %>%
    filter(species == species_name)
  
  # Calculate metrics
  metrics <- list(
    species = species_name,
    n_records = nrow(species_data),
    # Extent of occurrence (bounding box area)
    # eoo = st_area(st_convex_hull(st_union(species_data))),
    # Centroid of distribution
    # centroid = st_centroid(st_union(species_data)),
    # Range size (various methods available)
    # Nearest neighbour distances
    timestamp = Sys.time()
  )
  
  return(metrics)
}

#' Calculate all spatial metrics for all species
#'
calculate_all_metrics <- function(occurrences) {
  
  species_list <- unique(occurrences$species)
  message(glue("Calculating spatial metrics for {length(species_list)} species"))
  
  # metrics <- map(species_list, ~calculate_range_metrics(occurrences, .x))
  # metrics_df <- bind_rows(metrics)
  
  # return(metrics_df)
  
  return(NULL)
}

# ==============================================================================
# Density Analysis
# ==============================================================================

#' Create kernel density surface
#'
#' @param occurrences sf point object
#' @param bandwidth Kernel bandwidth in map units
#' @param resolution Grid resolution for output raster
#' @return SpatRaster object
#'
create_density_surface <- function(occurrences, bandwidth, resolution) {
  
  message("Creating kernel density surface...")
  
  # NOTE: Implement density estimation
  # This can be done using spatstat, terra, or other spatial packages
  # Example conceptual workflow:
  # 1. Convert to point pattern object
  # 2. Calculate kernel density
  # 3. Convert to raster/SpatRaster
  
  # density_raster <- ...
  
  # return(density_raster)
  
  return(NULL)
}

# ==============================================================================
# Spatial Visualisation
# ==============================================================================

#' Create map of occurrence points
#'
#' @param occurrences sf object
#' @param title Plot title
#' @return tmap object
#'
map_occurrences <- function(occurrences, title = "Species Occurrences") {
  
  message("Creating occurrence map...")
  
  # Create map using tmap
  # map <- tm_shape(occurrences) +
  #   tm_dots(col = "species", size = 0.1, alpha = 0.6) +
  #   tm_layout(title = title) +
  #   tm_scale_bar() +
  #   tm_compass()
  
  # return(map)
  
  return(NULL)
}

# ==============================================================================
# Execute Spatial Processing
# ==============================================================================

message("\n=== Starting Spatial Processing ===\n")

# Convert to spatial format
# occurrences_sf <- convert_to_sf(clean_data)

# Project if needed for distance-based analyses
# if (study_parameters$crs_projected != study_parameters$crs_wgs84) {
#   occurrences_projected <- project_spatial_data(
#     occurrences_sf, 
#     study_parameters$crs_projected
#   )
# }

# Create spatial grid and aggregate
# grid <- create_spatial_grid(
#   bbox = study_parameters$bbox,
#   cell_size = 10000,  # 10km cells (adjust as needed)
#   crs = study_parameters$crs_projected
# )
# 
# aggregated <- aggregate_to_grid(occurrences_projected, grid)

# Calculate spatial metrics
# spatial_metrics <- calculate_all_metrics(occurrences_sf)

# Create density surfaces (optional)
# density_surface <- create_density_surface(
#   occurrences_projected,
#   bandwidth = 50000,  # 50km bandwidth (adjust as needed)
#   resolution = 1000   # 1km resolution
# )

# ==============================================================================
# Save Spatial Outputs
# ==============================================================================

# Save as GeoJSON (GEE-compatible)
# if (options$save_geojson) {
#   geojson_file <- here(paths$processed, "occurrences_spatial.geojson")
#   st_write(occurrences_sf, geojson_file, delete_dsn = TRUE)
#   message(glue("GeoJSON saved to: {geojson_file}"))
# }

# Save aggregated grid
# grid_file <- here(paths$processed, "spatial_grid_aggregated.geojson")
# st_write(aggregated, grid_file, delete_dsn = TRUE)

# Save spatial metrics
# metrics_file <- here(paths$processed, create_filename("spatial_metrics", "rds"))
# saveRDS(spatial_metrics, metrics_file)

# Save density raster
# if (!is.null(density_surface)) {
#   raster_file <- here(paths$processed, "density_surface.tif")
#   writeRaster(density_surface, raster_file, overwrite = TRUE)
# }

# ==============================================================================
# Create Summary Maps
# ==============================================================================

# Generate maps for visual inspection
# occurrence_map <- map_occurrences(occurrences_sf)

# Save map
# map_file <- here(paths$figures, "occurrence_map.png")
# tmap_save(occurrence_map, map_file, width = 10, height = 8, units = "in", dpi = 300)
# message(glue("Map saved to: {map_file}"))

message("\n=== Script Notes ===")
message("This is an annotated placeholder script.")
message("Implement spatial processing using sf, terra, and related packages.")
message("Choose appropriate grid sizes, bandwidths, and metrics for your study.")
message("Consider computational efficiency for large datasets.")

# ==============================================================================
# End of Script
# ==============================================================================
