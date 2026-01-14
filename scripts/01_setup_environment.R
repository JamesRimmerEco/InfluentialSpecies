# ==============================================================================
# Script: 01_setup_environment.R
# Purpose: Install and load required R packages for the Influential Species 
#          Mapping project
# Author: [Author name]
# Date: [Date]
# ==============================================================================

# Description:
# This script sets up the R environment by installing and loading all required
# packages. Run this script once at the beginning of the project or when 
# setting up on a new machine.

# ==============================================================================
# Package Installation
# ==============================================================================

# Define required packages
required_packages <- c(
  # Data manipulation and tidyverse
  "tidyverse",      # Suite of data science packages
  "dplyr",          # Data manipulation
  "readr",          # Reading data files
  "tidyr",          # Data tidying
  
  # Spatial data handling
  "sf",             # Simple features for vector data
  "terra",          # Raster and vector spatial data
  "raster",         # Legacy raster support (if needed)
  
  # Biodiversity data
  "rgbif",          # Access to GBIF occurrence data
  "rinat",          # iNaturalist data access (optional)
  
  # Data validation and quality control
  "CoordinateCleaner", # Cleaning occurrence data
  "scrubr",         # Data cleaning utilities
  
  # Utilities
  "here",           # Project-relative file paths
  "glue",           # String interpolation
  "janitor",        # Data cleaning functions
  
  # Visualisation
  "ggplot2",        # Plotting
  "tmap"            # Thematic maps
)

# Function to install missing packages
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    message("Installing missing packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages, dependencies = TRUE)
  } else {
    message("All required packages are already installed.")
  }
}

# Install missing packages
install_if_missing(required_packages)

# ==============================================================================
# Load Libraries
# ==============================================================================

message("Loading required libraries...")

# Suppress package startup messages for cleaner output
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(rgbif)
  library(CoordinateCleaner)
  library(here)
  library(glue)
  library(janitor)
})

# ==============================================================================
# Session Information
# ==============================================================================

# Document package versions for reproducibility
session_info <- sessionInfo()
message("R version: ", session_info$R.version$version.string)
message("Platform: ", session_info$platform)

# Save session information to file
writeLines(
  capture.output(sessionInfo()),
  here("docs", "session_info.txt")
)

message("Environment setup complete!")
message("Session information saved to docs/session_info.txt")

# ==============================================================================
# Optional: renv for dependency management
# ==============================================================================

# Uncomment the following lines to use renv for package management:
# if (!requireNamespace("renv", quietly = TRUE)) {
#   install.packages("renv")
# }
# renv::init()    # Initialise renv for this project
# renv::snapshot() # Save current package versions

# ==============================================================================
# End of Script
# ==============================================================================
