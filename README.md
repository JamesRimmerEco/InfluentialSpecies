# Influential Species Mapping

## Project Overview

This repository supports exploratory research investigating influential species within a European context. It provides reproducible R-based workflows to:

- pull species occurrence records from GBIF and NBN,
- harmonise and merge per-species datasets across sources,
- apply conservative de-duplication,
- add QC flags for later filtering and diagnostics,
- prepare downstream products for modelling and rasterisation workflows (including Google Earth Engine handoff).

**Status**: Exploratory / short-term research project  
**Primary language**: R  
**Spatial platform**: Google Earth Engine integration (data preparation only)

## Objectives

1. Acquire and process species occurrence data from open data sources
2. Clean and check biodiversity records using transparent, staged rules
3. Prepare datasets for modelling / rasterisation workflows
4. Maintain reproducible and well-documented pipelines (regeneratable outputs)

## Repository Structure

```
InfluentialSpecies/
├── R/ # Core pipeline functions
├── scripts/ # Run scripts / wrappers (species sets, stage runners)
├── data/
│ ├── raw/ # Raw pull outputs (typically gitignored)
│ │ ├── gbif/<group>/<slug>/ # Per-species GBIF outputs (gbif_<slug>clean.csv)
│ │ └── nbn/<group>/<slug>/ # Per-species NBN outputs (nbn<slug>_clean.csv)
│ ├── processed/ # Stage outputs (typically gitignored; regeneratable)
│ │ ├── 01_merged/ # Stage 01 merged outputs
│ │ └── 02_qc_flagged/ # Stage 02 QC-flagged outputs
│ ├── outputs/ # Figures, tables, exports
│ └── credentials.R # Local GBIF creds (gitignored; optional; see below)
├── data/_checkpoints/ # Checkpoints for resumable pulls (gitignored)
│ ├── gbif/ # GBIF paging state + download keys (per species)
│ └── nbn/ # NBN cached raw pulls (per species)
├── docs/ # Documentation and notes
└── README.md # This file
```


## Pipeline

Stage outputs live under `data/processed/<stage>/...`. Outputs are designed to be regeneratable; large data files are typically gitignored, with optional `.gitkeep` placeholders.

### Stage 00 — Pull raw occurrences

- Engine: `R/pull_raw_occurrences.R`
- Wrappers: `scripts/pull_raw_species_set_*.R`

Outputs (per species):
- `data/raw/gbif/<group>/<slug>/gbif_<slug>_clean.csv`
- `data/raw/nbn/<group>/<slug>/nbn_<slug>_clean.csv`

Notes:
- NBN records should include a certainty grade where available.
- Uncertainty tags and related fields should be retained through ingestion so they can drive downstream filtering (especially geographic uncertainty and sensitive species handling). :contentReference[oaicite:1]{index=1}

### Stage 01 — Merge + conservative de-duplication

- Engine: `R/merge_dedup_occurrences.R`
- Wrapper example: `scripts/merge_dedup_species_set_6sp_test.R`

Outputs:
- `data/processed/01_merged/<slug>/occ_<slug>__merged.(parquet|rds)`
- Runlog: `data/processed/01_merged/_runlog_01_merged.csv`

Behaviour:
- Reads per-source CSVs (GBIF + NBN) with all columns forced to character to avoid type drift.
- Produces a merged per-species table with stable column retention.
- Performs conservative 1-to-1 cross-source duplicate removal (day + rounded coordinates), preferring a specified source where duplicates occur.
- Supports grouped and ungrouped Stage 01 layouts (controlled via `group_dir` in wrappers).

### Stage 02 — QC flagging (annotation only)

- Engine: `R/qc_flag_occurrences.R`
- Wrapper example: `scripts/qc_flag_species_set_6sp_test.R`

Outputs:
- `data/processed/02_qc_flagged/<slug>/occ_<slug>__qc_flagged.(parquet|rds)`
- Runlog: `data/processed/02_qc_flagged/_runlog_02_qc_flagged.csv`

Behaviour:
- Reads Stage 01 merged per-species outputs (parquet or rds).
- Adds QC flag columns to support later filtering and diagnostics.
- Does not drop records at this stage; it only annotates them.
- Stage 02 outputs are always written ungrouped under `data/processed/02_qc_flagged/<slug>/` (no test subfolders).

Flag set (current):
- `qc_flag_missing_coords`       : lon or lat is missing
- `qc_flag_coords_out_of_range`  : lon/lat present but outside valid ranges
- `qc_flag_missing_date`         : no day-level date and no year available
- `qc_flag_future_date`          : day-level date in the future, or year in the future
- `qc_flag_uncertainty_high`     : coordinateUncertaintyInMeters > threshold (default 10,000 m)
- `qc_flag_unexpected_licence`   : `licence_expected` indicates an unexpected licence (if present)
- `qc_flag_has_issues`           : `issues` field is non-empty (GBIF-style issues list)
- `qc_flag_any`                  : TRUE if any flag is TRUE
- `qc_flag_count`                : number of TRUE flags per record (optional helper)

## GBIF Record Limits and “Run Twice” Behaviour

GBIF’s standard search API is hard-limited to **100,000 records per query**. The pipeline handles this automatically:

- If a species has **≤ 100,000** GBIF records (for the configured query), the script uses a paged search and writes results immediately to `data/raw/gbif/...`.
- If a species has **> 100,000** GBIF records, the script switches to the GBIF **download** service (asynchronous):
  - The first run submits a download job and saves a **download key** to `data/_checkpoints/gbif/`.
  - Re-run the same wrapper script later to fetch and process the completed download.
  - The wrapper prints a final warning listing any species still pending.

## Caching / Re-running

Once a species is complete (the per-species CSV exists and its checkpoint is marked complete), re-running a Stage 00 wrapper will **not re-download** that species. Cached outputs are reused unless inputs are deleted or cache settings are changed.

Stage 01 and Stage 02 wrappers can be configured to skip cached outputs unless inputs are newer (see wrapper settings like `overwrite` and `refresh_if_inputs_newer`).

## Credentials (GBIF downloads only)

GBIF downloads require credentials. This repo supports a local, gitignored credentials file at:

- `data/credentials.R`

This file should define `GBIF_USER`, `GBIF_PWD`, and `GBIF_EMAIL`. Alternatively, set these as environment variables.

## Data Management

- Raw data (`data/raw/`) should never be modified directly.
- Processed data (`data/processed/`) contains stage outputs and should be treated as regeneratable.
- Outputs (`data/outputs/`) stores results, figures, and exports.
- Large data files should be stored externally when needed.

## Later stages

Further filtering, thinning/rasterisation, and modelling integration are planned downstream after Stage 02.
