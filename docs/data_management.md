# Data Management Guidelines

## Overview

This document outlines data management practices for the Influential Species Mapping project, including data organisation, storage, versioning, and sharing protocols.

## Data Directory Structure

```
data/
├── raw/                    # Original, unmodified data
│   ├── occurrences_raw_*.rds
│   └── [other raw downloads]
├── processed/              # Cleaned and processed data
│   ├── occurrences_clean_*.rds
│   ├── spatial_grid_*.geojson
│   └── project_config.rds
└── outputs/                # Analysis outputs and exports
    ├── figures/            # Plots and maps
    ├── tables/             # Summary tables
    └── gee_exports/        # GEE-ready files
```

## Data Categories

### Raw Data
**Location**: `data/raw/`

**Description**: Original, unmodified data as downloaded from source

**Principles**:
- **Never modify** raw data files directly
- Treat as read-only archives
- Document data provenance and download dates
- Keep separate from processed data

**Examples**:
- GBIF occurrence downloads
- Reference datasets
- External environmental layers

**Version Control**: Generally **not** tracked in Git (use `.gitignore`)

**Storage**: 
- Local: Project directory
- Backup: External drive or institutional storage
- Large files: Consider external repositories (e.g., Zenodo, Figshare)

---

### Processed Data
**Location**: `data/processed/`

**Description**: Cleaned, validated, and analysis-ready datasets

**Principles**:
- Fully reproducible from raw data + scripts
- Document processing steps
- Include quality flags and metadata
- Version using timestamps or sequential numbers

**Examples**:
- Cleaned occurrence records
- Validated spatial features
- Aggregated summaries
- Configuration files

**Version Control**: Small files can be tracked; use Git LFS for larger files

**Naming Convention**: `description_[timestamp].extension`

---

### Outputs
**Location**: `data/outputs/`

**Description**: Results, visualisations, and export-ready files

**Principles**:
- Reproducible from processed data + scripts
- Timestamped for version tracking
- Organised by type (figures, tables, exports)

**Examples**:
- Maps and plots
- Summary statistics tables
- GEE export files
- Manuscripts figures

**Version Control**: Generally not tracked (regenerated as needed)

## File Naming Conventions

### General Principles
- Use lowercase with underscores
- Be descriptive but concise
- Include dates/timestamps for versioning
- Use standard file extensions

### Examples

**Good**:
- `occurrences_clean_20240115.rds`
- `spatial_grid_10km_v2.geojson`
- `species_richness_map.png`

**Avoid**:
- `data.csv` (not descriptive)
- `FinalFINAL_v3_REAL.rds` (unclear versioning)
- `my analysis-output (1).shp` (spaces, special characters)

### Timestamps
Use ISO 8601 format: `YYYYMMDD` or `YYYYMMDD_HHMMSS`

Example: `occurrences_raw_20240115_143022.rds`

## Data Formats

### Recommended Formats

| Data Type | Format | Rationale |
|-----------|--------|-----------|
| Tabular data | CSV, RDS | Widely compatible; RDS preserves R objects |
| Spatial vector | GeoJSON, GeoPackage | Open standards; GEE compatible |
| Spatial raster | GeoTIFF | Standard format; retains metadata |
| Configuration | JSON, YAML | Human-readable; version-controllable |
| Documentation | Markdown, plain text | Accessible; version-controllable |

### Avoid
- Proprietary formats (e.g., .xlsx, .docx) for archival
- Compressed formats without documentation
- Binary formats without open specifications

## Metadata

### Essential Metadata
Each dataset should include:
- **Description**: What the data represents
- **Date created/modified**: When was it generated
- **Source**: Where did it come from
- **Processing**: How was it created
- **Contact**: Who to ask about the data
- **Licence**: Terms of use

### Implementation

**Option 1**: Separate metadata files
- Create `.json` or `.txt` metadata file
- Same name as data file with `_metadata` suffix
- Example: `occurrences_clean.rds` + `occurrences_clean_metadata.json`

**Option 2**: Embedded metadata
- Store in R object attributes
- Include in README files
- Document in notebooks/scripts

### Example Metadata File

```json
{
  "filename": "occurrences_clean_20240115.rds",
  "description": "Cleaned species occurrence records",
  "date_created": "2024-01-15T14:30:00Z",
  "source": "GBIF (www.gbif.org)",
  "processing_script": "04_clean_validate_data.R",
  "n_records": 45231,
  "n_species": 12,
  "crs": "EPSG:4326",
  "quality_checks": ["coordinate_validation", "duplicate_removal"],
  "contact": "researcher@institution.ac.uk"
}
```

## Storage and Backup

### Local Storage
- Primary location: Project directory on research computer
- Regular backups to external drive
- Sync to institutional cloud storage

### External Storage
**When to use**:
- Files > 100 MB
- Shared datasets
- Long-term archival
- Published data

**Options**:
- Institutional repository
- Zenodo (https://zenodo.org/)
- Figshare (https://figshare.com/)
- Dryad (https://datadryad.org/)
- Google Earth Engine Assets

### Version Control with Git

**Track in Git**:
- All scripts
- Documentation
- Small configuration files (< 1 MB)
- Metadata files

**Exclude from Git** (add to `.gitignore`):
- Raw data files
- Large processed files (> 50 MB)
- Temporary files
- Outputs (unless essential)

### Git LFS
For larger files that need version control:
```bash
git lfs track "*.rds"
git lfs track "*.tif"
```

## Data Size Management

### Strategies for Large Datasets

**Spatial Data**:
- Simplify geometries if appropriate
- Reduce coordinate precision (6 decimal places sufficient for most ecology)
- Aggregate to coarser resolution
- Use spatial indexing

**Tabular Data**:
- Remove unnecessary columns
- Use efficient data types (integer vs numeric)
- Compress with `gzip` for archival
- Split large files by species or region

**Raster Data**:
- Use appropriate compression (e.g., LZW)
- Reduce bit depth if possible
- Tile large rasters
- Consider Cloud Optimised GeoTIFF (COG)

## Data Sharing

### Internal Sharing
- Use institutional cloud storage (OneDrive, Google Drive)
- Share processed data, not raw (ensure reproducibility)
- Include README with instructions
- Version shared files clearly

### External Sharing
- Deposit in recognised repository
- Obtain DOI for citation
- Include comprehensive metadata
- Specify data licence (e.g., CC0, CC-BY)
- Follow FAIR principles (Findable, Accessible, Interoperable, Reusable)

### Data Licences
Consider:
- **CC0**: Public domain, no restrictions
- **CC-BY**: Attribution required
- **CC-BY-SA**: Attribution + share-alike

## Data Provenance

### Documentation
Track the origin and processing of all data:

```
occurrences_clean_20240115.rds
    ← 04_clean_validate_data.R
    ← occurrences_raw_20240115.rds
        ← 03_acquire_occurrence_data.R
        ← GBIF download DOI: 10.15468/dl.xxxxx
```

### Tools
- Git commits document script changes
- Timestamps track file versions
- Metadata files record processing steps
- Session info captures package versions

## Security and Privacy

### Sensitive Data
This project uses open biodiversity data, but consider:
- **Species of conservation concern**: May require location obfuscation
- **Personally identifiable information**: Remove if present
- **Proprietary data**: Ensure proper permissions and licences

### Best Practices
- Review data licences before sharing
- Check for sensitive species (IUCN Red List)
- Coordinate rounding for rare/threatened species
- Consult with data providers on sharing restrictions

## Quality Assurance

### Data Validation Checklist
- [ ] File names follow convention
- [ ] Metadata documented
- [ ] Processing steps recorded
- [ ] Quality checks performed
- [ ] Backups created
- [ ] Git commits informative
- [ ] File sizes appropriate
- [ ] Provenance clear

### Regular Maintenance
- Weekly: Check backups
- Monthly: Review directory organisation
- End of project: Deposit in repository, archive outdated files

## Tools and Resources

### R Packages
- `here`: Project-relative paths
- `fs`: File system operations
- `jsonlite`: JSON metadata
- `tools`: File utilities

### External Tools
- Git LFS: Large file versioning
- rsync: Backup synchronisation
- tar/gzip: File compression

### Resources
- GBIF Data Quality Requirements: https://www.gbif.org/data-quality-requirements
- FAIR Principles: https://www.go-fair.org/fair-principles/
- Data Carpentry: https://datacarpentry.org/

## Questions?

For data management queries, consult:
1. This documentation
2. Institutional data management guidelines
3. Project supervisor or data manager
4. Research data management support services

---

**Document Version**: 1.0  
**Last Updated**: [To be updated]  
**Maintained by**: Project team
