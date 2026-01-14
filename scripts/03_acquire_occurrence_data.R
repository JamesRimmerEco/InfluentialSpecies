# ==============================================================================
# Script: 03_acquire_occurrence_data.R
# Purpose: Download species occurrence records from biodiversity databases
# Author: [Author name]
# Date: [Date]
# ==============================================================================

# Description:
# This script acquires species occurrence data from open biodiversity databases,
# primarily using GBIF (Global Biodiversity Information Facility). The script
# handles API queries, data downloads, and initial formatting. Adjust the data
# source and query parameters based on your research needs.

# ==============================================================================
# Setup
# ==============================================================================

# Load required libraries
library(tidyverse)
library(rgbif)
library(here)
library(glue)

# Source project configuration
source(here("scripts", "02_configure_project.R"))

# ==============================================================================
# GBIF Data Acquisition Functions
# ==============================================================================

#' Download occurrence data for a single species from GBIF
#'
#' @param species_name Character string of scientific name
#' @param bbox List with xmin, xmax, ymin, ymax coordinates
#' @param date_range Vector with start and end dates
#' @return Dataframe of occurrence records
#' 
acquire_species_occurrences <- function(species_name, 
                                       bbox, 
                                       date_range) {
  
  message(glue("Querying GBIF for: {species_name}"))
  
  # NOTE: This is a placeholder function structure
  # Implement actual GBIF query using rgbif package functions such as:
  # - occ_search() for smaller datasets
  # - occ_download() for larger datasets (requires GBIF account)
  
  # Example structure (not functional without proper implementation):
  # occurrences <- occ_search(
  #   scientificName = species_name,
  #   geometry = paste(bbox$xmin, bbox$ymin, bbox$xmax, bbox$ymax, sep = ","),
  #   hasCoordinate = TRUE,
  #   hasGeospatialIssue = FALSE,
  #   limit = 10000
  # )
  
  message(glue("  Query complete for {species_name}"))
  
  # Return placeholder - replace with actual implementation
  return(NULL)
}

# ==============================================================================
# Batch Download for Multiple Species
# ==============================================================================

#' Download occurrence data for multiple species
#'
#' @param species_vector Character vector of species names
#' @return List of dataframes, one per species
#'
download_all_species <- function(species_vector) {
  
  message(glue("Beginning download for {length(species_vector)} species"))
  
  # Initialise storage
  all_occurrences <- list()
  
  # Loop through species
  for (i in seq_along(species_vector)) {
    
    species <- species_vector[i]
    message(glue("\n[{i}/{length(species_vector)}] Processing: {species}"))
    
    # Download data (replace with actual implementation)
    species_data <- acquire_species_occurrences(
      species_name = species,
      bbox = study_parameters$bbox,
      date_range = c(study_parameters$date_start, study_parameters$date_end)
    )
    
    # Store results
    all_occurrences[[species]] <- species_data
    
    # Polite delay between API calls (recommended by GBIF)
    Sys.sleep(1)
  }
  
  return(all_occurrences)
}

# ==============================================================================
# Execute Data Acquisition
# ==============================================================================

message("Starting species occurrence data acquisition...")
message(glue("Target species: {length(species_list)}"))
message(glue("Study extent: {study_parameters$bbox$xmin}, {study_parameters$bbox$ymin} to {study_parameters$bbox$xmax}, {study_parameters$bbox$ymax}"))
message(glue("Date range: {study_parameters$date_start} to {study_parameters$date_end}"))

# Download occurrence data
# NOTE: Uncomment and modify when ready to implement
# raw_occurrences <- download_all_species(species_list)

# ==============================================================================
# Format and Combine Data
# ==============================================================================

# Combine species data into single dataframe
# NOTE: Implement data formatting based on the structure returned by GBIF
# This typically includes:
# - Extracting relevant fields (species, coordinates, date, etc.)
# - Standardising column names
# - Converting data types
# - Adding metadata

# Example structure (not functional):
# combined_occurrences <- raw_occurrences %>%
#   bind_rows(.id = "species_query") %>%
#   select(
#     species = species,
#     longitude = decimalLongitude,
#     latitude = decimalLatitude,
#     date = eventDate,
#     basis_of_record = basisOfRecord,
#     coordinate_uncertainty = coordinateUncertaintyInMeters,
#     record_id = gbifID
#   ) %>%
#   mutate(
#     date = as.Date(date),
#     download_date = Sys.Date()
#   )

# ==============================================================================
# Save Raw Data
# ==============================================================================

# Save raw occurrence data
# output_file <- here(paths$raw_data, create_filename("occurrences_raw", "rds"))
# saveRDS(combined_occurrences, output_file)
# message(glue("Raw data saved to: {output_file}"))

# Also save as CSV for inspection
# if (options$save_csv) {
#   csv_file <- here(paths$raw_data, create_filename("occurrences_raw", "csv"))
#   write_csv(combined_occurrences, csv_file)
#   message(glue("CSV saved to: {csv_file}"))
# }

# ==============================================================================
# Summary Statistics
# ==============================================================================

# Generate and display summary statistics
# message("\n=== Data Acquisition Summary ===")
# message(glue("Total records downloaded: {nrow(combined_occurrences)}"))
# message(glue("Species with data: {length(unique(combined_occurrences$species))}"))
# message(glue("Date range: {min(combined_occurrences$date)} to {max(combined_occurrences$date)}"))

# Records per species summary
# species_summary <- combined_occurrences %>%
#   group_by(species) %>%
#   summarise(
#     n_records = n(),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(n_records))
# 
# print(species_summary)

message("\n=== Script Notes ===")
message("This is an annotated placeholder script.")
message("Implement GBIF API calls using rgbif package functions.")
message("Consider using occ_download() for large datasets (requires GBIF credentials).")
message("Ensure compliance with GBIF data use requirements and citation practices.")

# ==============================================================================
# End of Script
# ==============================================================================
