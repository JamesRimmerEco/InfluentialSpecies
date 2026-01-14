# Google Earth Engine Integration

## Overview

This document provides guidance on integrating R-processed species occurrence data with Google Earth Engine (GEE) for spatial analysis and visualisation. The workflow connects R-based data processing with GEE's cloud computing capabilities.

## Workflow Overview

```
R Processing → Export Files → Upload to GEE → GEE Analysis
(scripts 01-06)   (GeoJSON)    (Asset Manager)  (JavaScript)
```

## Prerequisites

### GEE Account
- Sign up at https://earthengine.google.com/signup/
- Academic/research accounts typically approved within 24-48 hours
- Requires Google account

### Required Tools
- **R/RStudio**: For data processing (scripts 01-06)
- **GEE Code Editor**: Web-based IDE (https://code.earthengine.google.com/)
- **earthengine-api** (optional): Python command-line tool

### Knowledge Requirements
- Basic JavaScript (for GEE scripting)
- R programming (covered in project scripts)
- GIS concepts (coordinate systems, vector/raster data)

## Data Preparation in R

### Output from R Pipeline
Script `06_prepare_gee_export.R` produces:
- **GeoJSON files**: Vector data in GEE-compatible format
- **Metadata files**: Documentation of data structure and content
- **Upload instructions**: Step-by-step guidance

### Data Format Requirements

**Coordinate System**: WGS84 (EPSG:4326)
- GEE uses WGS84 for all vector data
- R script handles conversion automatically

**File Size Limits**:
- Maximum 250 MB per asset
- For larger datasets: split by species, region, or time period

**Attribute Naming**:
- Use alphanumeric characters and underscores
- Avoid spaces and special characters
- Keep names concise (<20 characters preferred)

**Geometry Types**:
- Point: Individual occurrence records
- Polygon: Aggregated grid cells, study areas
- LineString: Transects, boundaries (if applicable)

## Uploading Assets to GEE

### Method 1: GEE Code Editor (Recommended for Beginners)

1. **Open Code Editor**: https://code.earthengine.google.com/

2. **Navigate to Assets Tab** (left panel)

3. **Upload New Asset**:
   - Click **"NEW"** button
   - Select **"Table Upload"** (for vector data)

4. **Select Files**:
   - Choose your `.geojson` file from `data/outputs/gee_exports/`
   - GEE supports: GeoJSON, Shapefiles (zipped), CSV with coordinates

5. **Configure Upload**:
   - **Asset ID**: Create meaningful name (e.g., `influential_species_points`)
   - **Properties**: Review attribute columns
   - **CRS**: Should auto-detect WGS84

6. **Start Upload**:
   - Click "Upload"
   - Monitor progress in Tasks tab

7. **Wait for Processing**:
   - Upload time varies with file size
   - Check Tasks tab for completion status

### Method 2: Python API (For Batch Uploads)

```python
import ee

# Initialise Earth Engine
ee.Initialize()

# Define upload parameters
asset_id = 'users/YOUR_USERNAME/influential_species_occurrences'
file_path = '/path/to/file.geojson'

# Create upload task
ee.batch.Export.table.toAsset(
    collection=ee.FeatureCollection(file_path),
    description='Upload influential species data',
    assetId=asset_id
).start()

print(f"Upload started: {asset_id}")
```

### Method 3: Command Line (earthengine-api)

```bash
# Install earthengine-api
pip install earthengine-api

# Authenticate
earthengine authenticate

# Upload asset
earthengine upload table \
  --asset_id=users/YOUR_USERNAME/influential_species_points \
  /path/to/occurrences.geojson

# Check status
earthengine task list
```

## Using Assets in GEE Scripts

### Accessing Uploaded Assets

```javascript
// Load your uploaded asset
var speciesOccurrences = ee.FeatureCollection('users/YOUR_USERNAME/influential_species_points');

// Inspect properties
print('Number of records:', speciesOccurrences.size());
print('First record:', speciesOccurrences.first());

// Get property names
var properties = speciesOccurrences.first().propertyNames();
print('Properties:', properties);
```

### Basic Visualisation

```javascript
// Define visualisation parameters
var visParams = {
  color: 'red',
  pointSize: 3,
  pointShape: 'circle',
  width: 1,
  fillColor: 'ff000033'  // Semi-transparent red
};

// Add to map
Map.addLayer(speciesOccurrences, visParams, 'Species Occurrences');
Map.centerObject(speciesOccurrences, 6);
```

### Filtering Data

```javascript
// Filter by species
var targetSpecies = speciesOccurrences.filter(
  ee.Filter.eq('species', 'Species name')
);

// Filter by date range
var recentRecords = speciesOccurrences.filter(
  ee.Filter.date('2020-01-01', '2023-12-31')
);

// Filter by region
var studyArea = ee.Geometry.Rectangle([-10, 49, 2, 61]);  // UK bbox
var regionalRecords = speciesOccurrences.filterBounds(studyArea);
```

### Spatial Analysis Examples

#### Species Richness Heatmap

```javascript
// Create grid cells
var grid = studyArea.coveringGrid(ee.Projection('EPSG:4326'), 10000);  // 10km

// Count species per cell
var speciesRichness = grid.map(function(cell) {
  var count = speciesOccurrences.filterBounds(cell).size();
  return cell.set('richness', count);
});

// Visualise
var richnessVis = {
  min: 0,
  max: 50,
  palette: ['blue', 'yellow', 'red']
};
Map.addLayer(speciesRichness, richnessVis, 'Species Richness');
```

#### Occurrence Density

```javascript
// Create density surface using reduceToImage
var density = speciesOccurrences
  .reduceToImage({
    properties: ['species'],
    reducer: ee.Reducer.countDistinct()
  })
  .reproject('EPSG:4326', null, 1000);  // 1km resolution

// Visualise
Map.addLayer(
  density, 
  {min: 0, max: 10, palette: ['white', 'blue', 'green', 'yellow', 'red']},
  'Density'
);
```

#### Temporal Animation

```javascript
// Group by year
var years = ee.List.sequence(2000, 2023);

var yearlyCollections = years.map(function(year) {
  var yearFilter = ee.Filter.calendarRange(year, year, 'year');
  return speciesOccurrences
    .filter(yearFilter)
    .set('year', year);
});

// Create animation (see GEE documentation for full implementation)
```

## Integrating with GEE Datasets

### Environmental Data

```javascript
// Load climate data
var temperature = ee.ImageCollection('ECMWF/ERA5/DAILY')
  .select('mean_2m_air_temperature')
  .filterDate('2020-01-01', '2020-12-31')
  .mean();

// Sample at occurrence points
var sampledData = temperature.sampleRegions({
  collection: speciesOccurrences,
  scale: 1000,
  geometries: true
});

print('Sampled data:', sampledData.limit(5));
```

### Land Cover Analysis

```javascript
// Load land cover
var landCover = ee.Image('ESA/WorldCover/v100/2020');

// Extract land cover at occurrences
var occurrencesWithLC = landCover.sampleRegions({
  collection: speciesOccurrences,
  properties: ['species'],
  scale: 10
});

// Summarise by land cover type
var lcSummary = occurrencesWithLC.aggregate_histogram('Map');
print('Land cover distribution:', lcSummary);
```

## Exporting Results from GEE

### Export to Google Drive

```javascript
// Export processed data
Export.table.toDrive({
  collection: processedData,
  description: 'species_analysis_results',
  fileFormat: 'GeoJSON'
});

// Export raster
Export.image.toDrive({
  image: densityRaster,
  description: 'density_map',
  scale: 1000,
  region: studyArea
});
```

### Export to Asset (for Reuse)

```javascript
Export.table.toAsset({
  collection: processedData,
  description: 'processed_occurrences',
  assetId: 'users/YOUR_USERNAME/processed_occurrences'
});
```

## Best Practices

### Asset Organisation
- Use consistent naming: `project_datatype_version`
- Create folders for different projects
- Document asset IDs in R project documentation
- Include date/version in asset descriptions

### Computational Efficiency
- **Filter early**: Reduce data before expensive operations
- **Use appropriate scale**: Match analysis resolution to data
- **Batch processing**: Use `map()` for repetitive operations
- **Cache results**: Export intermediate results as assets

### Data Quality
- Verify CRS after upload
- Check attribute preservation
- Test with small subset first
- Validate geometry integrity

### Code Documentation
```javascript
// Always document your GEE scripts
// Include:
// - Purpose of the script
// - Asset IDs used
// - Analysis parameters
// - Expected outputs
// - Date and author
```

## Troubleshooting

### Upload Issues

**Error**: "Invalid geometry"
- **Solution**: Check coordinate validity in R, ensure WGS84

**Error**: "File too large"
- **Solution**: Split data spatially or temporally, aggregate points to grid

**Error**: "Asset not found"
- **Solution**: Check asset path, ensure upload completed

### Processing Issues

**Error**: "Computation timed out"
- **Solution**: Reduce spatial extent, simplify operations, export intermediate results

**Error**: "Memory limit exceeded"
- **Solution**: Process in batches, use `.filterBounds()`, reduce resolution

## Resources

### GEE Documentation
- **Earth Engine Guide**: https://developers.google.com/earth-engine/guides
- **API Reference**: https://developers.google.com/earth-engine/apidocs
- **Tutorials**: https://developers.google.com/earth-engine/tutorials
- **Data Catalogue**: https://developers.google.com/earth-engine/datasets

### Community Resources
- **GEE Developers Forum**: https://groups.google.com/g/google-earth-engine-developers
- **Stack Overflow**: Tag `google-earth-engine`
- **GitHub Examples**: https://github.com/google/earthengine-api

### Learning Materials
- **GEE JavaScript Tutorial**: https://developers.google.com/earth-engine/tutorials/tutorial_js_01
- **Python API Guide**: https://developers.google.com/earth-engine/guides/python_install
- **Geospatial Course**: https://courses.spatialthoughts.com/

## Example GEE Script Template

```javascript
// ============================================================================
// Script: Influential Species GEE Analysis
// Purpose: Analyse species occurrences with environmental data
// Author: [Your name]
// Date: [Date]
// ============================================================================

// --- Configuration ---
var assetPath = 'users/YOUR_USERNAME/influential_species_occurrences';
var studyRegion = ee.Geometry.Rectangle([-10, 49, 2, 61]);
var dateStart = '2000-01-01';
var dateEnd = '2023-12-31';

// --- Load Data ---
var occurrences = ee.FeatureCollection(assetPath);

print('Total occurrences:', occurrences.size());
print('First record:', occurrences.first());

// --- Filter Data ---
var filtered = occurrences
  .filterBounds(studyRegion)
  .filter(ee.Filter.date(dateStart, dateEnd));

// --- Analysis ---
// [Add your analysis code here]

// --- Visualisation ---
Map.centerObject(filtered, 6);
Map.addLayer(filtered, {color: 'red'}, 'Occurrences');

// --- Export ---
// [Add export tasks if needed]
```

## Integration Workflow Summary

1. **Process data in R** (scripts 01-06)
2. **Export GEE-ready files** (script 06)
3. **Upload to GEE** (Code Editor or API)
4. **Develop GEE scripts** (JavaScript)
5. **Run analyses** (GEE cloud computing)
6. **Export results** (back to Drive or as assets)
7. **Further analysis in R** (if needed)

## Questions?

For GEE-specific questions:
- Check GEE documentation first
- Search GEE Developers Forum
- Post detailed question with reproducible example

For R-GEE integration issues:
- Review this document
- Check script 06 output
- Verify file formats and CRS

---

**Document Version**: 1.0  
**Last Updated**: [To be updated]  
**GEE Version**: Current (updated regularly by Google)
