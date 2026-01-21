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
├── R/                          # Core R functions (GBIF/NBN pulls, cleaning utilities)
├── scripts/                    # Run scripts / wrappers (e.g. species-set pulls)
├── data/
│   ├── raw/                    # Raw pull outputs (typically gitignored)
│   │   ├── gbif/<group>/<slug>/ # Per-species GBIF outputs (gbif_<slug>_clean.csv)
│   │   └── nbn/<group>/<slug>/  # Per-species NBN outputs (nbn_<slug>_clean.csv)
│   ├── processed/              # Cleaned / harmonised datasets for analysis
│   ├── outputs/                # Figures, tables, exports
│   └── credentials.R           # Local GBIF creds (gitignored; optional; see docs)
├── data/_checkpoints/          # Checkpoints for resumable pulls (gitignored)
│   ├── gbif/                   # GBIF paging state + download keys (per species)
│   └── nbn/                    # NBN cached raw pulls (per species)
├── docs/                       # Documentation and notes
└── README.md                   # This file
```

## Getting Started
### GBIF + NBN data pulls (important notes)

Species occurrence pulls are handled by `R/pull_raw_occurrences.R` and are typically run via a wrapper script in `scripts/` (e.g. `scripts/pull_raw_species_set_*.R`).

### GBIF record limits and “run twice” behaviour

GBIF’s standard search API is hard-limited to **100,000 records per query**. The pipeline handles this automatically:

- If a species has **≤ 100,000** GBIF records (for the configured query), the script uses a paged search and writes results immediately to `data/raw/gbif/...`.
- If a species has **> 100,000** GBIF records, the script switches to the GBIF **download** service (asynchronous):
  - The first run will submit a download job and save a **download key** to `data/_checkpoints/gbif/`.
  - You may need to **run the same wrapper script again later** to fetch and process the completed download.
  - The script prints a final **GBIF WARNING** listing any species still pending.

### Caching / re-running

Once a species is complete (the per-species CSV exists and its checkpoint is marked complete), re-running the wrapper script will **not re-download** that species. It will reuse cached outputs unless you delete the files or set `use_cache = FALSE`.

### Credentials (GBIF downloads only)

GBIF downloads require credentials. This repo supports a local, gitignored credentials file at:

- `data/credentials.R` (not synced)

This file should define `GBIF_USER`, `GBIF_PWD`, and `GBIF_EMAIL`. Alternatively, you can set these as environment variables.

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

