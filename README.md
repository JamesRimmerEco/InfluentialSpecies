# Influential Species Mapping

## Project Overview

This repository supports exploratory postdoctoral research investigating influential species within ecological networks. The project focuses on reproducible R-based workflows for processing species occurrence data, conducting spatial analyses, and preparing datasets for integration with Google Earth Engine (GEE).

**Status**: Exploratory / Short-term research project  
**Primary Language**: R  
**Spatial Platform**: Google Earth Engine integration

## Objectives

1. Acquire and process species occurrence data from open data sources
2. Clean and validate biodiversity records for spatial analysis
3. Conduct spatial processing and summary statistics
4. Prepare spatially-explicit datasets for Google Earth Engine workflows
5. Maintain reproducible and well-documented analytical pipelines

## Repository Structure

```
InfluentialSpecies/
├── R/                      # R functions and utilities
├── scripts/                # Analysis scripts (numbered workflow)
├── data/
│   ├── raw/               # Original, unmodified data (not tracked)
│   ├── processed/         # Cleaned and processed data
│   └── outputs/           # Analysis outputs and exports
├── docs/                  # Documentation and notes
└── README.md              # This file
```

## Getting Started

### Prerequisites

- R (≥ 4.0.0)
- RStudio (recommended)
- Required R packages (see `scripts/01_setup_environment.R`)

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

3. Configure data paths and parameters in `scripts/02_configure_project.R`

## Workflow

The analysis workflow is designed to be run sequentially:

1. **01_setup_environment.R** - Install and load required packages
2. **02_configure_project.R** - Set project parameters and paths
3. **03_acquire_occurrence_data.R** - Download species occurrence records
4. **04_clean_validate_data.R** - Clean and validate occurrence data
5. **05_spatial_processing.R** - Conduct spatial analyses and summaries
6. **06_prepare_gee_export.R** - Format data for Google Earth Engine

See `docs/workflow.md` for detailed documentation of each step.

## Data Management

- **Raw data** (`data/raw/`) should never be modified directly
- **Processed data** (`data/processed/`) contains cleaned, analysis-ready datasets
- **Outputs** (`data/outputs/`) stores results, figures, and GEE-ready exports
- Large data files should be stored externally (see `docs/data_management.md`)

## Google Earth Engine Integration

This project prepares data for use in Google Earth Engine but does not include GEE JavaScript code. See `docs/gee_integration.md` for guidance on:

- Exporting data in GEE-compatible formats
- Uploading assets to GEE
- Linking R outputs with GEE workflows

## Contributing

This is an exploratory research project. If you wish to contribute or collaborate:

1. Open an issue to discuss proposed changes
2. Fork the repository and create a feature branch
3. Ensure code is well-documented and follows existing style
4. Submit a pull request with a clear description

See `CONTRIBUTING.md` for detailed guidelines.

## Reproducibility

To ensure reproducibility:

- All scripts are version-controlled
- Package versions are documented
- Random seeds are set where applicable
- Data provenance is tracked in metadata files

## Citation

If you use this repository in your research, please cite:

```
[Citation information to be added upon publication]
```

## Licence

[Licence to be determined - consult with institution and collaborators]

## Contact

For questions or collaboration enquiries, please open an issue or contact the repository maintainer.

## Acknowledgements

This research was conducted as part of postdoctoral work exploring ecological network dynamics and species influence patterns.
