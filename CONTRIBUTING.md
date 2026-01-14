# Contributing to Influential Species Mapping

Thank you for your interest in contributing to this project! This document provides guidelines for contributing to the Influential Species Mapping research repository.

## Project Context

This is an exploratory postdoctoral research project focused on reproducible R-based workflows for species occurrence data analysis. Contributions should align with the project's academic research goals and maintain scientific rigour.

## How to Contribute

### Reporting Issues

If you encounter bugs, data quality concerns, or have suggestions:

1. **Search existing issues** to avoid duplicates
2. **Open a new issue** with a clear, descriptive title
3. **Provide details**:
   - What you expected to happen
   - What actually happened
   - Steps to reproduce (if applicable)
   - Your R version and package versions (use `sessionInfo()`)
   - Any error messages or warnings

### Suggesting Enhancements

For new features or improvements:

1. **Open an issue** to discuss the proposal
2. **Explain the rationale**: Why is this enhancement needed?
3. **Describe the approach**: How might it be implemented?
4. **Consider implications**: Impact on existing workflows?

### Contributing Code

#### Before You Start

1. **Discuss significant changes** via an issue first
2. **Check existing code style** and follow conventions
3. **Ensure changes are minimal** and focused
4. **Consider reproducibility** implications

#### Development Process

1. **Fork the repository**
2. **Create a feature branch** (`git checkout -b feature/your-feature-name`)
3. **Make your changes**:
   - Follow existing code style
   - Add comments where helpful
   - Update documentation if needed
4. **Test your changes**:
   - Run affected scripts
   - Check for errors or warnings
   - Verify reproducibility
5. **Commit your changes**:
   - Use clear, descriptive commit messages
   - Reference related issues
6. **Push to your fork**
7. **Submit a pull request**

#### Pull Request Guidelines

**PR Title**: Clear, concise description of changes

**PR Description** should include:
- **Purpose**: What problem does this solve?
- **Changes**: What did you modify?
- **Testing**: How did you verify the changes?
- **Documentation**: Did you update docs?
- **Breaking changes**: Any compatibility issues?

**Code Review**:
- Be open to feedback
- Respond to comments promptly
- Make requested changes in new commits

## Code Style

### R Style Guidelines

Follow the [tidyverse style guide](https://style.tidyverse.org/):

**Naming**:
- Functions: `lowercase_with_underscores`
- Variables: `descriptive_names`
- Constants: `UPPERCASE_WITH_UNDERSCORES`

**Spacing**:
```r
# Good
x <- 5
result <- function_name(arg1, arg2)

# Avoid
x<-5
result <- function_name( arg1,arg2 )
```

**Line Length**: Maximum 80 characters when possible

**Comments**:
- Explain *why*, not *what*
- Use `#` for inline comments
- Use `#'` for function documentation

**Code Structure**:
```r
# Section headers with ====
# ==============================================================================
# Section Name
# ==============================================================================

# Subsections with ----
# Subsection name ----

# Regular comments for explanation
```

### Documentation

- **Scripts**: Include header with purpose, author, date
- **Functions**: Document parameters and return values
- **Changes**: Update relevant documentation files
- **Examples**: Provide where helpful

### File Organisation

- Place R functions in `R/` directory
- Place analysis scripts in `scripts/` directory
- Number scripts sequentially (01_, 02_, etc.)
- Use descriptive file names

## Data Contributions

### Adding New Data Sources

If proposing new data sources:
1. Document data provenance
2. Check licencing and permissions
3. Follow data management guidelines
4. Update `docs/data_management.md`

### Data Quality

- Ensure data quality checks are implemented
- Document any data cleaning steps
- Maintain data provenance records

## Testing

While this project doesn't have formal unit tests, please:
- **Run affected scripts** after changes
- **Check output** for expected results
- **Verify reproducibility** on clean environment
- **Document test procedures** in PR

## Documentation

Update documentation when making changes:
- **README.md**: High-level project information
- **docs/workflow.md**: Workflow changes
- **docs/data_management.md**: Data handling changes
- **docs/gee_integration.md**: GEE-related changes
- **Script comments**: Inline documentation

## Versioning

This project uses:
- **Git tags** for major milestones
- **Timestamps** for data file versions
- **Commit messages** for change tracking

## Reproducibility Requirements

All contributions must maintain reproducibility:
- Use relative paths (via `here` package)
- Document package dependencies
- Set random seeds where applicable
- Avoid hard-coded absolute paths
- Record R and package versions

## Communication

### Preferred Channels

- **GitHub Issues**: Bug reports, feature requests
- **Pull Requests**: Code contributions
- **Discussions**: General questions (if enabled)

### Response Times

This is a research project with limited resources:
- Issues/PRs may take several days for response
- Complex changes may require extended discussion
- Patience and understanding are appreciated

## Code of Conduct

### Expected Behaviour

- Be respectful and inclusive
- Welcome newcomers
- Accept constructive criticism
- Focus on project goals
- Acknowledge contributions

### Unacceptable Behaviour

- Harassment or discrimination
- Trolling or insulting comments
- Personal attacks
- Unprofessional conduct

### Enforcement

Violations should be reported to project maintainers. Responses may include warnings, temporary blocks, or permanent bans.

## Licence

By contributing, you agree that your contributions will be licenced under the same terms as the project.

## Attribution

Contributors will be acknowledged in:
- Git commit history
- Future publications (for significant contributions)
- CONTRIBUTORS.md file (if created)

## Questions?

If you're unsure about contributing:
1. Review existing issues and PRs
2. Read project documentation
3. Open an issue to ask
4. Contact maintainers

## Getting Help

### Resources

- **Project Documentation**: `/docs` directory
- **R Documentation**: `?function_name` or `help(package_name)`
- **Tidyverse**: https://www.tidyverse.org/
- **R Spatial**: https://r-spatial.org/

### Common Issues

**Environment Setup**:
- Run `01_setup_environment.R` first
- Check package versions match session info
- Verify R version compatibility

**Data Issues**:
- Check file paths use `here::here()`
- Verify data directory structure
- Review data management guidelines

**Git Issues**:
- Ensure `.gitignore` is respected
- Don't commit large data files
- Use meaningful commit messages

## Recognition

We appreciate all contributions, including:
- Code improvements
- Bug reports
- Documentation enhancements
- Workflow suggestions
- Data quality checks
- Script testing

Thank you for helping improve this research project!

---

**Document Version**: 1.0  
**Last Updated**: [To be updated]
