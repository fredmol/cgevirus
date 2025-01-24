# Changelog

## [1.3.0] - 2024-01-24
### Removed
- Removed CDD (Conserved Domain Detection) analysis and reporting functionality
- Removed CDD-related database dependencies
- Simplified report layout by removing CDD sections

## [1.2.0] - 2024-01-24
### Added
- PDF report generation with detailed analysis results
- Integration with WeasyPrint for PDF creation
- Custom DTU styling for reports
- GC content calculation for virus sequences

### Dependencies Added
- weasyprint>=63.0
- Jinja2>=3.1.4
- matplotlib>=3.9.2