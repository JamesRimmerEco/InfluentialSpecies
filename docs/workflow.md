# Workflow Documentation

## Overview

This document describes the analytical workflow for the Influential Species Mapping project. The workflow is designed as a sequential pipeline, with each script building upon the outputs of previous steps.

## Workflow Diagram

```
01_setup_environment.R
         ↓
02_configure_project.R
         ↓
03_acquire_occurrence_data.R → raw data
         ↓
04_clean_validate_data.R → processed data
         ↓
05_spatial_processing.R → spatial outputs
         ↓
06_prepare_gee_export.R → GEE-ready files
```

## Detailed Steps

### Step 1: Environment Setup
**Script**: `01_setup_environment.R`

**Purpose**: Install and load required R packages

**Inputs**: None

**Outputs**: 
- Installed R packages
- Session information file (`docs/session_info.txt`)

**Key Actions**:
- Checks for and installs missing packages
- Loads required libraries
- Documents package versions for reproducibility

**Run Once**: At project initiation or when setting up on a new machine

---

### Step 2: Project Configuration
**Script**: `02_configure_project.R`

**Purpose**: Define project parameters and settings

**Inputs**: User-defined parameters (within script)

**Outputs**: 
- Project configuration object (`data/processed/project_config.rds`)
- Created directory structure

**Key Actions**:
- Sets study area boundaries
- Defines species list
- Establishes quality control thresholds
- Creates output directories

**Modify**: Adjust parameters based on your specific research question

---

### Step 3: Occurrence Data Acquisition
**Script**: `03_acquire_occurrence_data.R`

**Purpose**: Download species occurrence records from GBIF

**Inputs**: 
- Species list (from Step 2)
- Study parameters (from Step 2)

**Outputs**: 
- Raw occurrence data (`data/raw/occurrences_raw_[timestamp].rds`)
- CSV backup (optional)

**Key Actions**:
- Queries GBIF API for each species
- Downloads occurrence records within study extent
- Filters by date range
- Saves raw data for archival

**Data Source**: Global Biodiversity Information Facility (GBIF)

**Note**: Requires internet connection. Large downloads may take considerable time.

---

### Step 4: Data Cleaning and Validation
**Script**: `04_clean_validate_data.R`

**Purpose**: Quality control and validation of occurrence records

**Inputs**: 
- Raw occurrence data (from Step 3)
- Quality parameters (from Step 2)

**Outputs**: 
- Cleaned occurrence data (`data/processed/occurrences_clean_[timestamp].rds`)
- Flagged records (`data/processed/occurrences_flagged_[timestamp].rds`)
- CSV backups (optional)

**Key Actions**:
- Removes records with missing coordinates
- Eliminates duplicates
- Validates coordinate accuracy
- Checks temporal consistency
- Applies CoordinateCleaner tests
- Filters by basis of record
- Flags high uncertainty records

**Quality Checks**:
- Coordinate validity (range, format)
- Spatial outliers
- Country/sea mismatches
- Institutional/centroid coordinates
- Temporal validity
- Coordinate uncertainty thresholds

---

### Step 5: Spatial Processing
**Script**: `05_spatial_processing.R`

**Purpose**: Spatial analysis and aggregation

**Inputs**: 
- Cleaned occurrence data (from Step 4)

**Outputs**: 
- Spatial features (`data/processed/occurrences_spatial.geojson`)
- Aggregated grid (`data/processed/spatial_grid_aggregated.geojson`)
- Spatial metrics (`data/processed/spatial_metrics_[timestamp].rds`)
- Density surfaces (`data/processed/density_surface.tif`)
- Maps (`data/outputs/figures/occurrence_map.png`)

**Key Actions**:
- Converts to sf spatial objects
- Projects to appropriate CRS
- Creates spatial grids
- Aggregates points to cells
- Calculates spatial metrics (range size, density)
- Generates kernel density surfaces
- Produces summary maps

**Spatial Analyses**:
- Point-to-grid aggregation
- Species richness by grid cell
- Extent of occurrence (EOO)
- Range centroids
- Kernel density estimation

---

### Step 6: GEE Export Preparation
**Script**: `06_prepare_gee_export.R`

**Purpose**: Format data for Google Earth Engine upload

**Inputs**: 
- Spatial data (from Step 5)

**Outputs**: 
- GEE-formatted files (`data/outputs/gee_exports/*.geojson`)
- Asset metadata (`data/outputs/gee_exports/*_metadata.json`)
- Upload instructions (`data/outputs/gee_exports/GEE_UPLOAD_INSTRUCTIONS.md`)

**Key Actions**:
- Rounds coordinates to specified precision
- Simplifies attribute tables
- Converts dates to GEE format
- Optimises file sizes
- Generates metadata
- Creates upload documentation

**GEE Requirements**:
- WGS84 coordinate system (EPSG:4326)
- Maximum file size: 250 MB per asset
- Compatible attribute names (no special characters)
- Appropriate geometry types (Point, Polygon)

---

## Running the Workflow

### Complete Workflow
To run the entire workflow from start to finish:

```r
# Step 1: Setup (run once)
source("scripts/01_setup_environment.R")

# Step 2: Configure (modify parameters first)
source("scripts/02_configure_project.R")

# Step 3: Acquire data
source("scripts/03_acquire_occurrence_data.R")

# Step 4: Clean data
source("scripts/04_clean_validate_data.R")

# Step 5: Spatial processing
source("scripts/05_spatial_processing.R")

# Step 6: Prepare for GEE
source("scripts/06_prepare_gee_export.R")
```

### Partial Workflow
To run from a specific step (e.g., after modifying cleaning parameters):

```r
# Load configuration
source("scripts/02_configure_project.R")

# Re-run cleaning with new parameters
source("scripts/04_clean_validate_data.R")

# Continue downstream steps
source("scripts/05_spatial_processing.R")
source("scripts/06_prepare_gee_export.R")
```

## Customisation Points

### Study Area
Modify `study_parameters$bbox` in `02_configure_project.R`

### Target Species
Edit `species_list` in `02_configure_project.R`

### Quality Thresholds
Adjust `quality_parameters` in `02_configure_project.R`

### Spatial Resolution
Change grid cell size in `05_spatial_processing.R`

### GEE Format
Set `options$gee_format` in `02_configure_project.R`

## Error Handling

Common issues and solutions:

**API Rate Limits**: Add delays between GBIF queries or use download requests for large datasets

**Memory Issues**: Process species individually or use data subsetting

**Coordinate System Errors**: Ensure consistent CRS across all spatial operations

**Large File Sizes**: Aggregate spatially or split into multiple assets

## Data Flow

```
Raw Data (GBIF) 
    → Raw occurrences (unmodified archive)
    → Cleaned occurrences (quality controlled)
    → Spatial features (georeferenced)
    → Aggregated summaries (gridded)
    → GEE assets (formatted for upload)
```

## Quality Assurance

At each step, review:
- Console messages for warnings/errors
- Summary statistics
- Flagged records
- Output file sizes
- Visual plots

## Version Control

- All scripts are version-controlled via Git
- Track parameter changes in commit messages
- Use timestamped filenames for outputs
- Document major analytical decisions

## Reproducibility Checklist

- [ ] Session information documented
- [ ] Package versions recorded
- [ ] Random seeds set
- [ ] Parameters logged in config file
- [ ] Data provenance tracked
- [ ] Metadata files generated
- [ ] Scripts commented

## Next Steps After This Workflow

1. Upload assets to Google Earth Engine
2. Develop GEE JavaScript analysis code
3. Integrate R outputs with GEE workflows
4. Conduct ecological analyses
5. Generate publication-ready figures
6. Document findings

## References

- GBIF.org: https://www.gbif.org/
- CoordinateCleaner: https://github.com/ropensci/CoordinateCleaner
- Google Earth Engine: https://earthengine.google.com/
- sf package: https://r-spatial.github.io/sf/

## Contact

For questions about this workflow, open an issue in the repository or contact the project maintainer.
