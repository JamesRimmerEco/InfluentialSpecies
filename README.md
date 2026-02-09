# InfluentialSpecies

A reproducible (and deliberately transparent) pipeline for pulling, harmonising, QC-flagging, policy-filtering, and gridding species occurrence records for downstream analysis (e.g. SDMs / MaxEnt).

This README reflects the **cloud_safe** branch. It is set up so you can run long pulls on a machine that syncs the repo via Google Drive/OneDrive without committing credentials or relying on fragile checkpoints stored inside the repo.

---

## What this repo does

The pipeline is organised into stages:

- Stage 00 — Pull “raw clean” occurrences from GBIF (Europe-wide) + NBN Atlas (UK), with only:
  - coordinate screening (must have lon/lat)
  - light-touch within-source de-duplication (exact duplicates only)

- Stage 01 — Merge GBIF+NBN into a single table per species and remove obvious 1-to-1 duplicates across sources.

- Stage 02 — Add QC flags (no hard filtering yet): uncertainty, missingness, obvious problems, etc.

- Stage 03 — Apply an explicit policy filter (this is where “what counts as acceptable evidence” is defined).

- Stage 04 — Grid/rasterise Stage 03 occurrences (e.g. 25 km cells), with optional land masking and quick map/context outputs.

Intent:
- Stages 00–02 do not “throw records away” beyond minimal sanity checks.
- Stage 03 is the first time we make policy decisions about what is admissible.

---

## Repository structure

Only code + lightweight docs are tracked. Raw and processed data are written under `data/` and ignored by git.

Pipeline-style layout (high level):

├── R/                        # Core pipeline functions (stage “engines”, helpers)
│
├── scripts/                  # Run scripts / wrappers (species sets, stage runners)
│
├── data/
│
│  ├── raw/                   # Raw pull outputs (typically gitignored)
│  │
│  │  ├── gbif/<group>/<slug>/  # Per-species GBIF outputs (e.g. gbif_<slug>_clean.csv)
│  │  └── nbn/<group>/<slug>/   # Per-species NBN outputs  (e.g. nbn_<slug>_clean.csv)
│  │
│  ├── processed/             # Stage outputs (typically gitignored; regeneratable)
│  │
│  │  ├── 01_merged/          # Stage 01 merged outputs (GBIF + NBN combined; light de-dupe)
│  │  ├── 02_qc_flagged/      # Stage 02 QC-flagged outputs (flags added; little/no dropping)
│  │  ├── 03_policy_filtered/ # Stage 03 policy-filtered outputs (records retained per policy)
│  │  ├── 04_grid/            # Stage 04 gridded outputs (e.g. 25 km EPSG:3035)
│  │  │  └── <grid_run_tag>/<slug>/   # Per-run tag (extent/landmask/resolution) then species
│  │  │     ├── occ_<slug>__grid25km.parquet
│  │  │     ├── presence_centroids_25km_<slug>.csv    # lon,lat,cell_id,n_points_in_cell
│  │  │     └── ... (optional per-species diagnostics)
│  │  └── 05_maxent/          # Stage 05 SDM/MaxEnt outputs (optional local outputs/exports)
│  │     └── <model_run_tag>/<slug>/  # e.g. predictors/resolution/region tags
│  │
│  ├── outputs/               # Figures, tables, exports (e.g. summary CSVs, maps)
│  │
│  └── credentials.R          # Local GBIF creds (gitignored; optional; see docs)
│
├── data/_checkpoints/        # Checkpoints for resumable pulls (gitignored)
│
│  ├── gbif/                  # GBIF paging state + download keys (per species)
│  └── nbn/                   # NBN cached raw pulls (per species)
│
├── docs/                     # Documentation and notes
│
└── README.md                 # Project overview

Notes:
- `<group_dir>` is a run folder used for Stage 00 pulls (e.g. `home_run_2026-02-06`).
- `<slug>` is a filesystem-safe species name (lowercase, underscores).
- `<policy_id>` / `<policy_tag>` are identifiers for Stage 03/04 configuration.

---

## Prerequisites

### R packages

Core (Stage 00):
- rgbif, galah, dplyr, stringr, readr, lubridate, tibble

Stages 01–04 additionally use:
- arrow (Parquet I/O), data.table
- Stage 04 map output additionally: sf, rnaturalearth, rnaturalearthdata, ggplot2

Install (example):

  install.packages(c(
    "rgbif","galah","dplyr","stringr","readr","lubridate","tibble",
    "arrow","data.table","sf","rnaturalearth","rnaturalearthdata","ggplot2"
  ))

---

## Credentials and logins (cloud_safe)

### GBIF (needed for >100k species downloads)

The code looks for credentials in the environment first:
- GBIF_USER
- GBIF_PWD
- GBIF_EMAIL

Recommended: put these in your user-level `.Renviron` (not in the repo).

There is also a local credentials-file fallback that is deliberately gitignored:
- credentials.R (repo root), or
- data/credentials.R

If credentials are already present in the environment, the engine will not override them.

### NBN Atlas (galah)

Typical login:

  library(galah)
  galah_config(atlas = "United Kingdom", email = "you@domain", verbose = FALSE)
  galah_login()

NBN tokens can expire mid-run; wrappers attempt “best effort” login but you may need to re-run `galah_login()` interactively if you start seeing HTTP 403 / OAuth errors.

---

## Disk and temp-file hygiene (important on Windows)

GBIF downloads for common birds/mammals can be hundreds of MB to multiple GB per species (zips + extracted occurrence tables). On Windows this can silently fill `C:` because:
- checkpoints may live under LOCALAPPDATA (depending on wrapper settings)
- R’s `tempdir()` lives under `C:\Users\<you>\AppData\Local\Temp`
- if R aborts mid-unzip/read, temporary folders can be left behind

Two controls are recommended.

### 1) Put checkpoints (including GBIF download zips) on a drive with space

Set this before running wrappers:

  Sys.setenv(INFLUENTIAL_CHECKPOINT_ROOT = "G:/InfluentialSpecies_checkpoints")
  dir.create(Sys.getenv("INFLUENTIAL_CHECKPOINT_ROOT"), recursive = TRUE, showWarnings = FALSE)

This moves:
- GBIF download zips to .../gbif/
- NBN checkpoints to .../nbn/

### 2) Move R temporary files off C: (recommended for huge downloads)

Best practice is to set Windows environment variables (System Properties → Environment Variables) so that TEMP and TMP point to a large drive (e.g. `G:\Temp`), then restart R/RStudio.

If you want a session-level attempt:

  Sys.setenv(TMPDIR = "G:/TempR", TEMP = "G:/TempR", TMP = "G:/TempR")
  dir.create("G:/TempR", recursive = TRUE, showWarnings = FALSE)

Note: R often chooses `tempdir()` at session start, so the most reliable method is setting Windows env vars and restarting.

---

## Pipeline transparency (policy choices)

See:
- docs/pipeline_transparency.qmd

This document explains:
- what each stage does and does not do
- where policy decisions occur (Stage 03)
- how policy settings are recorded in outputs and runlogs

If you use Quarto:

  quarto render docs/pipeline_transparency.qmd

---

## How to run the pipeline

### Stage 00 — Pull raw occurrences (GBIF + NBN)

Engine:
- R/pull_raw_occurrences.R

Wrapper templates:
- scripts/pull_raw_species_set_6sp_test_v_2.0.R
  A small test run. Useful to check credentials, paths, and that the engines work end-to-end.

- scripts/pull_raw_species_set_mapping_list_100sp_HOME_SAFE_v1.R
  A long “home safe” run intended for unattended execution over a large species list.
  Designed to resume cleanly in a single stable group_dir.

Key behaviour (GBIF):
- GBIF occ_search is capped at 100,000 records per query.
- The engine checks expected record count:
  - <= 100k: uses occ_search paging and writes the CSV immediately
  - > 100k: uses GBIF downloads (occ_download), which are asynchronous
- For >100k species, the first run typically:
  - submits a download job
  - saves the download key to a checkpoint
  - continues to the next species
- Re-running later resumes downloads when they are ready and writes the final CSV.

Outputs (Stage 00):
- data/raw/gbif/<group_dir>/gbif_<slug>_clean.csv
- data/raw/nbn/<group_dir>/nbn_<slug>_clean.csv

---

### Stage 01 — Merge GBIF + NBN and remove obvious cross-source duplicates

Engine:
- R/merge_dedup_occurrences.R

Wrapper example:
- scripts/merge_dedup_species_set_6sp_test.R

Outputs:
- data/processed/01_merged/<slug>/occ_<slug>__merged.parquet
- data/processed/01_merged/_runlog_01_merged.csv

This stage is conservative:
- merges sources and removes only obvious duplicates
- does not apply policy filters

---

### Stage 02 — QC flagging (no hard filtering)

Engine:
- R/qc_flag_occurrences.R

Wrapper example:
- scripts/qc_flag_species_set_6sp_test.R

Outputs:
- data/processed/02_qc_flagged/<slug>/occ_<slug>__stage02_qc_flagged.parquet
- data/processed/02_qc_flagged/_runlog_02_qc_flagged.csv

Adds QC flags used for Stage 03 filtering and diagnostics.

---

### Stage 03 — Policy filtering (explicit choices)

Engine:
- R/filter_occurrences.R

Wrapper example:
- scripts/filter_species_set_6sp_test.R

Stage 03 requires a named policy list (documented in docs/pipeline_transparency.qmd).
Policies typically specify (examples):
- year window
- maximum coordinate uncertainty
- accepted record types (e.g. in situ observations)
- source restrictions (GBIF-only / NBN-only / both)
- licence handling rules (if required)

Outputs:
- data/processed/03_filtered/<policy_id>/<slug>/occ_<slug>__stage03_filtered.parquet
- data/processed/03_filtered/_runlog_03_filtered.csv

---

### Stage 04 — Gridding / rasterisation of Stage 03 outputs

Engine:
- R/grid_occurrences_stage03.R (main function: grid_stage03_to_grid())

Wrapper example:
- scripts/grid_stage03_25km_test.R

Outputs:
- data/processed/04_grid/<policy_tag>/_summary_grid.csv
- data/processed/04_grid/<policy_tag>/<slug>_presence_cells.csv
- optional: data/processed/04_grid/<policy_tag>/_bbox_context_map.png

---

## Example “new user” workflow

1) Clone repo and switch to the branch:

  git clone <repo-url>
  cd InfluentialSpecies
  git checkout cloud_safe

2) In R, set a checkpoint root on a large drive (recommended):

  Sys.setenv(INFLUENTIAL_CHECKPOINT_ROOT = "G:/InfluentialSpecies_checkpoints")
  dir.create(Sys.getenv("INFLUENTIAL_CHECKPOINT_ROOT"), recursive = TRUE, showWarnings = FALSE)

3) Run a small Stage 00 test:

  source("scripts/pull_raw_species_set_6sp_test_v_2.0.R")

4) Run Stage 01–04 using the 6-species wrappers as templates:

  source("scripts/merge_dedup_species_set_6sp_test.R")
  source("scripts/qc_flag_species_set_6sp_test.R")
  source("scripts/filter_species_set_6sp_test.R")
  source("scripts/grid_stage03_25km_test.R")

5) For long runs (e.g. 100 species), use:

  source("scripts/pull_raw_species_set_mapping_list_100sp_HOME_SAFE_v1.R")

That wrapper is designed to:
- run unattended for long periods
- resume in a single stable group_dir
- skip species that already have both outputs in that group_dir
- tolerate transient failures (per-species tryCatch)
- resume GBIF downloads on later runs

---

## Troubleshooting (common)

### “It’s re-pulling a species I already have”

Stage 00 will only fully skip a species when:
- the clean CSV exists, and
- the checkpoint is marked complete = TRUE

If R crashed mid-run, you can end up with:
- a CSV on disk but complete = FALSE in the checkpoint

In that case the engine may “verify completeness” and re-pull. If you are confident the CSV is complete, you can mark the checkpoint complete (or re-run once and let it repair itself).

### NBN 403 / OAuth errors mid-run

Re-authenticate:

  library(galah)
  galah_config(atlas = "United Kingdom", email = "you@domain", verbose = FALSE)
  galah_login()

### “My C: drive keeps filling up”

Use the disk hygiene controls above:
- set INFLUENTIAL_CHECKPOINT_ROOT to a large drive
- move TEMP/TMP off C: (restart R afterwards)

If you need to reclaim space:
- delete large GBIF zip files under your checkpoint root: <checkpoint_root>/gbif/*.zip
- clear abandoned temp unzip folders under your R temp directory (often under LOCALAPPDATA\Temp)

---

## Notes

- This repo intentionally separates:
  - pull + harmonise (Stages 00–02)
  - policy decisions (Stage 03)
  - gridding/rasterisation configuration (Stage 04)

- If you change a Stage 03 policy definition, treat it as a new configuration and write outputs to a new policy_id / policy_tag.
