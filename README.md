# Influential Species Mapping

## Project Overview

This repository supports exploratory research investigating influential species within a European context. This project area focuses on reproducible R-based workflows for accessing and processing species occurrence data,  and preparing datasets for integration with Google Earth Engine (GEE).

**Status**: Exploratory / Short-term research project  
**Primary Language**: R  
**Spatial Platform**: Google Earth Engine integration

## Objectives

1. Acquire and process species occurrence data from open data sources
2. Clean and check biodiversity records
3. Prepare datasets for Google Earth Engine workflows
4. Maintain reproducible and well-documented analytical pipelines

## Repository Structure

```
InfluentialSpecies/
├── R/                      # R functions and utilities
├── scripts/                # Analysis scripts
├── data/
│   ├── raw/               # Original, unmodified data (not tracked)
│   ├── processed/         # Cleaned and processed data
│   └── outputs/           # Outputs and exports
├── docs/                  # Documentation and notes
└── README.md              # This file
```

## Getting Started

### Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/JamesRimmerEco/InfluentialSpecies.git
   cd InfluentialSpecies
   ```

2. Install required R packages:
   ```r
   source("scripts/01_setup_environment.R")
   ```

## Data Management

- **Raw data** (`data/raw/`) should never be modified directly
- **Processed data** (`data/processed/`) contains cleaned, analysis-ready datasets
- **Outputs** (`data/outputs/`) stores results, figures, and GEE-ready exports
- Large data files should be stored externally

## Google Earth Engine Integration

This project prepares data for use in Google Earth Engine but does not include GEE JavaScript code. 

## Reproducibility

To ensure reproducibility:

- All scripts are version-controlled
- Data provenance is tracked in metadata files

