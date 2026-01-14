# ==============================================================================
# Script: 04_clean_validate_data.R
# Purpose: Clean and validate species occurrence records
# Author: [Author name]
# Date: [Date]
# ==============================================================================

# Description:
# This script performs data quality control on species occurrence records,
# including coordinate validation, spatial outlier detection, duplicate removal,
# and taxonomic checking. The cleaning steps follow best practices for
# biodiversity data quality assessment.

# ==============================================================================
# Setup
# ==============================================================================

# Load required libraries
library(tidyverse)
library(sf)
library(CoordinateCleaner)
library(here)
library(glue)

# Source project configuration
source(here("scripts", "02_configure_project.R"))

# ==============================================================================
# Load Data
# ==============================================================================

message("Loading raw occurrence data...")

# Load the raw occurrence data from previous script
# raw_data <- readRDS(here(paths$raw_data, "occurrences_raw_[timestamp].rds"))

# For development/testing, create placeholder structure
# raw_data <- tibble(
#   species = character(),
#   longitude = numeric(),
#   latitude = numeric(),
#   date = Date(),
#   coordinate_uncertainty = numeric(),
#   basis_of_record = character(),
#   record_id = character()
# )

# ==============================================================================
# Data Cleaning Functions
# ==============================================================================

#' Remove records with missing coordinates
#'
remove_missing_coords <- function(data) {
  initial_n <- nrow(data)
  
  cleaned <- data %>%
    filter(!is.na(longitude), !is.na(latitude))
  
  removed <- initial_n - nrow(cleaned)
  message(glue("Removed {removed} records with missing coordinates"))
  
  return(cleaned)
}

#' Remove duplicate records
#'
remove_duplicates <- function(data) {
  initial_n <- nrow(data)
  
  cleaned <- data %>%
    distinct(species, longitude, latitude, date, .keep_all = TRUE)
  
  removed <- initial_n - nrow(cleaned)
  message(glue("Removed {removed} duplicate records"))
  
  return(cleaned)
}

#' Flag records with high coordinate uncertainty
#'
flag_uncertain_coords <- function(data, max_uncertainty) {
  data <- data %>%
    mutate(
      high_uncertainty = coordinate_uncertainty > max_uncertainty | 
                        is.na(coordinate_uncertainty)
    )
  
  n_flagged <- sum(data$high_uncertainty, na.rm = TRUE)
  message(glue("Flagged {n_flagged} records with uncertainty > {max_uncertainty}m"))
  
  return(data)
}

# ==============================================================================
# Coordinate Validation
# ==============================================================================

#' Perform coordinate validation using CoordinateCleaner
#'
#' This function applies multiple spatial checks including:
#' - Validity (coordinates within valid ranges)
#' - Country centroids
#' - Sea coordinates (for terrestrial species)
#' - Institutional coordinates
#' - Zero coordinates
#' - Urban areas (optional)
#'
validate_coordinates <- function(data) {
  
  message("Performing coordinate validation checks...")
  
  # NOTE: Implement CoordinateCleaner functions
  # Example workflow (requires actual implementation):
  
  # Convert to spatial object for CoordinateCleaner
  # data_sf <- st_as_sf(
  #   data, 
  #   coords = c("longitude", "latitude"),
  #   crs = study_parameters$crs_wgs84
  # )
  
  # Apply CoordinateCleaner tests
  # flags <- clean_coordinates(
  #   x = data,
  #   lon = "longitude",
  #   lat = "latitude",
  #   species = "species",
  #   tests = c("capitals", "centroids", "equal", "gbif", "institutions",
  #             "outliers", "seas", "zeros")
  # )
  
  # Add flags to original data
  # data_flagged <- data %>%
  #   mutate(
  #     passed_coord_checks = flags$.summary,
  #     flag_capitals = flags$.cap,
  #     flag_centroids = flags$.cen,
  #     flag_equal = flags$.equ,
  #     flag_gbif = flags$.gbf,
  #     flag_institutions = flags$.inst,
  #     flag_outliers = flags$.otl,
  #     flag_seas = flags$.sea,
  #     flag_zeros = flags$.zer
  #   )
  
  # Summary of flagged records
  # message(glue("Records passing all checks: {sum(data_flagged$passed_coord_checks)}"))
  # message(glue("Records flagged: {sum(!data_flagged$passed_coord_checks)}"))
  
  return(data)
}

# ==============================================================================
# Temporal Validation
# ==============================================================================

#' Validate and filter temporal data
#'
validate_temporal <- function(data, date_start, date_end) {
  
  message("Validating temporal data...")
  
  initial_n <- nrow(data)
  
  cleaned <- data %>%
    filter(
      !is.na(date),
      date >= as.Date(date_start),
      date <= as.Date(date_end),
      date <= Sys.Date()  # No future dates
    )
  
  removed <- initial_n - nrow(cleaned)
  message(glue("Removed {removed} records with invalid dates"))
  
  return(cleaned)
}

# ==============================================================================
# Basis of Record Filtering
# ==============================================================================

#' Filter by basis of record
#'
filter_basis_of_record <- function(data, acceptable_basis) {
  
  initial_n <- nrow(data)
  
  cleaned <- data %>%
    filter(basis_of_record %in% acceptable_basis)
  
  removed <- initial_n - nrow(cleaned)
  message(glue("Removed {removed} records with unacceptable basis of record"))
  
  return(cleaned)
}

# ==============================================================================
# Execute Cleaning Pipeline
# ==============================================================================

message("\n=== Starting Data Cleaning Pipeline ===\n")

# Apply cleaning steps sequentially
# cleaned_data <- raw_data %>%
#   remove_missing_coords() %>%
#   remove_duplicates() %>%
#   flag_uncertain_coords(quality_parameters$max_coord_uncertainty) %>%
#   validate_temporal(
#     study_parameters$date_start,
#     study_parameters$date_end
#   ) %>%
#   filter_basis_of_record(quality_parameters$acceptable_basis) %>%
#   validate_coordinates()

# ==============================================================================
# Filter Based on Quality Flags
# ==============================================================================

# Create final clean dataset by filtering flagged records
# final_clean <- cleaned_data %>%
#   filter(
#     passed_coord_checks == TRUE,
#     high_uncertainty == FALSE
#   )

# Retain flagged dataset for inspection
# flagged_records <- cleaned_data %>%
#   filter(
#     passed_coord_checks == FALSE |
#     high_uncertainty == TRUE
#   )

# ==============================================================================
# Summary Statistics
# ==============================================================================

# message("\n=== Cleaning Summary ===")
# message(glue("Initial records: {nrow(raw_data)}"))
# message(glue("Final clean records: {nrow(final_clean)}"))
# message(glue("Records flagged: {nrow(flagged_records)}"))
# message(glue("Retention rate: {round(nrow(final_clean)/nrow(raw_data)*100, 1)}%"))

# Per-species summary
# species_summary <- final_clean %>%
#   group_by(species) %>%
#   summarise(
#     n_records = n(),
#     date_range_start = min(date),
#     date_range_end = max(date),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(n_records))
# 
# print(species_summary)

# ==============================================================================
# Save Cleaned Data
# ==============================================================================

# Save final clean dataset
# output_file <- here(paths$processed, create_filename("occurrences_clean", "rds"))
# saveRDS(final_clean, output_file)
# message(glue("\nCleaned data saved to: {output_file}"))

# Save flagged records for inspection
# flagged_file <- here(paths$processed, create_filename("occurrences_flagged", "rds"))
# saveRDS(flagged_records, flagged_file)
# message(glue("Flagged records saved to: {flagged_file}"))

# Save as CSV if requested
# if (options$save_csv) {
#   csv_file <- here(paths$processed, create_filename("occurrences_clean", "csv"))
#   write_csv(final_clean, csv_file)
# }

message("\n=== Script Notes ===")
message("This is an annotated placeholder script.")
message("Implement cleaning functions using CoordinateCleaner and custom validation.")
message("Review flagged records to understand data quality issues.")
message("Adjust quality parameters based on your research requirements.")

# ==============================================================================
# End of Script
# ==============================================================================
