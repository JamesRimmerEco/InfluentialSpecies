# Getting Started

## Quick Start Guide

This guide will help you set up and run the Influential Species Mapping project.

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/JamesRimmerEco/InfluentialSpecies.git
cd InfluentialSpecies
```

### 2. Install R and RStudio

**R** (version ≥ 4.0.0):
- Download from: https://cran.r-project.org/
- Follow installation instructions for your OS

**RStudio** (recommended):
- Download from: https://posit.co/download/rstudio-desktop/
- Install after R is installed

### 3. Install Required Packages

Open RStudio and run:

```r
source("scripts/01_setup_environment.R")
```

This will install all required packages. First-time setup may take 10-20 minutes.

## Configuration

### 4. Configure Project Parameters

Edit `scripts/02_configure_project.R`:

```r
# Define your study area
study_parameters <- list(
  bbox = list(
    xmin = -10.0,  # Western longitude
    xmax = 2.0,    # Eastern longitude
    ymin = 49.0,   # Southern latitude
    ymax = 61.0    # Northern latitude
  ),
  date_start = "2000-01-01",
  date_end = "2023-12-31"
)

# Define your target species
species_list <- c(
  "Passer domesticus",  # Example: House Sparrow
  "Turdus merula",      # Example: Common Blackbird
  "Bombus terrestris"   # Example: Buff-tailed Bumblebee
)
```

Run the configuration script:

```r
source("scripts/02_configure_project.R")
```

## Running the Analysis

### Sequential Workflow

Run scripts in order:

```r
# 1. Setup (already done above)
source("scripts/01_setup_environment.R")

# 2. Configuration (already done above)
source("scripts/02_configure_project.R")

# 3. Acquire occurrence data
source("scripts/03_acquire_occurrence_data.R")

# 4. Clean and validate data
source("scripts/04_clean_validate_data.R")

# 5. Spatial processing
source("scripts/05_spatial_processing.R")

# 6. Prepare for Google Earth Engine
source("scripts/06_prepare_gee_export.R")
```

### RStudio Workflow

**Option 1**: Run line-by-line
- Open script in RStudio
- Place cursor on line
- Press `Ctrl+Enter` (Windows/Linux) or `Cmd+Enter` (Mac)

**Option 2**: Run entire script
- Open script in RStudio
- Click "Source" button or press `Ctrl+Shift+S`

**Option 3**: Run selection
- Highlight code section
- Press `Ctrl+Enter`

## Understanding Outputs

### Data Directories

After running the workflow, you'll have:

```
data/
├── raw/
│   └── occurrences_raw_[timestamp].rds    # Downloaded GBIF data
├── processed/
│   ├── occurrences_clean_[timestamp].rds  # Cleaned data
│   ├── occurrences_spatial.geojson        # Spatial format
│   └── spatial_metrics_[timestamp].rds    # Analysis results
└── outputs/
    ├── figures/
    │   └── occurrence_map.png             # Visualisations
    └── gee_exports/
        └── *.geojson                      # GEE-ready files
```

### Viewing Results

**Load processed data**:
```r
library(tidyverse)
library(here)

# Load cleaned occurrences
clean_data <- readRDS(here("data", "processed", "occurrences_clean_[timestamp].rds"))

# View summary
glimpse(clean_data)
summary(clean_data)
```

**View spatial data**:
```r
library(sf)
library(tmap)

# Load spatial data
spatial_data <- st_read(here("data", "processed", "occurrences_spatial.geojson"))

# Quick map
tm_shape(spatial_data) + tm_dots()
```

## Troubleshooting

### Common Issues

#### Package Installation Fails

**Issue**: Error installing packages

**Solutions**:
```r
# Update R to latest version
# Then try installing packages individually:
install.packages("tidyverse")
install.packages("sf")

# On Linux, you may need system dependencies:
# sudo apt-get install libudunits2-dev libgdal-dev libgeos-dev libproj-dev
```

#### GBIF Download Fails

**Issue**: Cannot download occurrence data

**Solutions**:
1. Check internet connection
2. Verify GBIF website is accessible: https://www.gbif.org/
3. Reduce number of species in species_list
4. Check API rate limits (add delays between queries)

#### Memory Errors

**Issue**: "Cannot allocate vector of size..."

**Solutions**:
```r
# Increase memory limit (Windows)
memory.limit(size = 16000)  # 16GB

# Process fewer species at once
# Work with data subsets
# Close other applications
```

#### Spatial Processing Slow

**Issue**: Scripts taking too long

**Solutions**:
- Start with small study area
- Reduce number of species
- Use coarser spatial resolution
- Process species individually

#### File Not Found Errors

**Issue**: "cannot open file..."

**Solutions**:
```r
# Always use here() for paths
library(here)
file_path <- here("data", "raw", "myfile.rds")

# Check working directory
getwd()

# Check file exists
file.exists(file_path)
```

## Customisation

### Modify Workflow for Your Needs

**Different study region**:
- Edit `study_parameters$bbox` in `02_configure_project.R`

**Different time period**:
- Edit `date_start` and `date_end` in `02_configure_project.R`

**Different data sources**:
- Modify `03_acquire_occurrence_data.R` to use other APIs (iNaturalist, etc.)

**Different quality thresholds**:
- Edit `quality_parameters` in `02_configure_project.R`

**Custom analyses**:
- Add new scripts following numbering system (e.g., `07_custom_analysis.R`)
- Add custom functions to `R/` directory

## Working with Scripts

### Script Structure

All scripts follow this structure:

```r
# ==============================================================================
# Script: [Name]
# Purpose: [Description]
# ==============================================================================

# Load libraries
library(package)

# Source configuration
source(here("scripts", "02_configure_project.R"))

# Load data
data <- readRDS(...)

# Analysis
results <- analyze(data)

# Save outputs
saveRDS(results, ...)

# ==============================================================================
# End of Script
# ==============================================================================
```

### Adding Your Own Code

1. **Small modifications**: Edit existing scripts
2. **New functions**: Add to `R/helper_functions.R` or create new file in `R/`
3. **New analyses**: Create new numbered script (e.g., `07_my_analysis.R`)

## Best Practices

### Data Management
- Never modify files in `data/raw/`
- Use timestamps in file names
- Document data processing steps
- Back up important results

### Code Style
- Use meaningful variable names
- Add comments for complex operations
- Follow tidyverse style guide
- Keep scripts modular and focused

### Reproducibility
- Always run `02_configure_project.R` first
- Use `set.seed()` for random operations
- Document R and package versions
- Keep scripts self-contained

## Next Steps

Once you've run the basic workflow:

1. **Explore data**: View cleaned occurrence data, check quality
2. **Visualise**: Create maps and plots
3. **Analyse**: Calculate spatial metrics, test hypotheses
4. **GEE Integration**: Upload to Google Earth Engine for further analysis
5. **Document**: Record findings in `docs/project_notes.md`

## Learning Resources

### R Programming
- **R for Data Science**: https://r4ds.had.co.nz/
- **Advanced R**: https://adv-r.hadley.nz/
- **RStudio Cheatsheets**: https://www.rstudio.com/resources/cheatsheets/

### Spatial Analysis in R
- **Geocomputation with R**: https://geocompr.robinlovelace.net/
- **Spatial Data Science**: https://r-spatial.github.io/sf/
- **rspatial**: https://rspatial.org/

### Biodiversity Data
- **GBIF**: https://www.gbif.org/
- **CoordinateCleaner**: https://ropensci.github.io/CoordinateCleaner/
- **rgbif**: https://docs.ropensci.org/rgbif/

### Google Earth Engine
- **GEE Guide**: https://developers.google.com/earth-engine/
- **GEE Tutorials**: https://developers.google.com/earth-engine/tutorials/

## Getting Help

1. **Check documentation**: Review files in `docs/` directory
2. **Search online**: R error messages are usually well-documented
3. **R help**: Use `?function_name` for function documentation
4. **Community**: Stack Overflow, RStudio Community
5. **Issues**: Open an issue in this repository

## Project Documentation

Complete documentation available in:
- `README.md`: Project overview
- `docs/workflow.md`: Detailed workflow documentation
- `docs/data_management.md`: Data handling guidelines
- `docs/gee_integration.md`: Google Earth Engine integration
- `CONTRIBUTING.md`: Contribution guidelines

## Feedback

If you encounter issues or have suggestions for improving this guide, please open an issue in the repository.

---

**Happy Analysing!**
