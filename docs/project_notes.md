# Project Notes

## Research Context

This repository supports exploratory postdoctoral research on influential species in ecological networks. The work is short-term and focuses on developing reproducible workflows rather than finalising specific ecological hypotheses.

## Current Status

**Phase**: Repository setup and workflow development  
**Last Updated**: [Date]

## Key Research Questions (Preliminary)

These are exploratory questions subject to refinement:

1. How are influential species distributed spatially?
2. What environmental factors correlate with influential species occurrence?
3. Can we identify spatial hotspots of influential species diversity?

## Methodological Considerations

### Data Sources
- **Primary**: GBIF (Global Biodiversity Information Facility)
- **Supplementary**: Environmental layers from Google Earth Engine
- **Temporal scope**: 2000-present (adjustable)

### Spatial Scope
- **Initial focus**: UK and surrounding regions (example)
- **Rationale**: Data availability, manageable scale for exploratory work
- **Flexibility**: Framework designed for any geographic region

### Species Selection
- Species list to be determined based on ecological network analysis
- Criteria for "influential species" still under development
- Placeholder species list in configuration script

### Analytical Approach
- **Descriptive**: Spatial patterns and distributions
- **Correlative**: Environmental associations
- **Reproducible**: All steps documented and scriptable

## Known Limitations

1. **Data Quality**: GBIF data has known biases (sampling effort, taxonomy)
2. **Temporal Coverage**: Variable across species and regions
3. **Spatial Resolution**: Limited by occurrence record precision
4. **Causation**: Analyses are correlative, not causal

## Future Directions

Potential extensions of this work:
- Integration with ecological network analysis
- Temporal trend analysis
- Species distribution modelling
- Climate change projections
- Conservation prioritisation

## Important Reminders

- This is exploratory work - methods may change
- Always validate assumptions
- Document decisions and rationale
- Keep raw data unmodified
- Maintain reproducibility at all stages

## Collaboration Notes

- Repository designed for solo postdoc work but open to collaboration
- Contact before major contributions
- See CONTRIBUTING.md for guidelines

## Technical Notes

### Performance Considerations
- Large GBIF downloads can take hours
- Spatial processing can be memory-intensive
- Consider working with data subsets during development

### Common Issues
- GBIF API rate limiting: Add delays between requests
- Memory errors: Process species individually or use sampling
- CRS mismatches: Always check and transform explicitly

## References and Resources

### Key Papers
[To be added as literature review progresses]

### Useful Resources
- GBIF Data Quality documentation
- CoordinateCleaner vignettes
- Google Earth Engine guides
- R Spatial resources

## Meeting Notes

[Space for meeting notes, decisions, feedback]

---

**Note**: This is a living document. Update regularly as the project evolves.
