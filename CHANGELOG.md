# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.2] UNRELEASED

### Bugfix

- Fixed broken hyperlinks in substructure/blast search results

## [1.1.1] 09-12-2024

### Added

- Implemented an API (/api/<mite_accession>) to serve MITE JSON files

### Changed

- Changed generic MITE contributor "AAAAA..." with explanation
- Improved search functionality for substructure
- Implemented generating of enzyme visualization from scratch

## [1.1.0] 08-12-2024

### Added

- Added search functionality (substructure and BLAST) to repository overview
- Added page for canonicalization or SMILES
- Added resolution of organism or origin
- Improved Download page

## [1.0.0] 30-11-24

### Added

- Added additional "tailoring" qualifiers

### Changed

- Fixed broken MIBiG links

## [0.2.8] 26-11-24

### Added

- Added page to allow generating SMILES strings from single-letter amino acid string
- Added a dark theme for the webpage

### Bugfix

- Fixed the setting of the 'isIntermediate' flag in the submission forms

## [0.2.7] 22-11-24

### Changed

- Updated Contact page
- Updated Help

## [0.2.6] 18-11-24

### Added

- Added information Troubleshooting
- Added cross-link to help/troubleshooting when error message is flashed
- Added video tutorial link
- Added jsondiff to show changes in modified files

## [0.2.5] 12-11-24

### Bugfix

- Fixed issue list link
- Added enzyme name to Issue title

## [0.2.4] 12-11-24

### Bugfix

- Changed the default reviewer ID to "BBBBB..." (before: default contributor ID "AAA...")

## [0.2.3] 12-11-24

### Added

- Added Help page
- Added script to upload submitted data as Issues to GitHub

### Changed

- Changed Dockerfile to automatically populate app with data (before: `prepare_mite_data.py` had to be run manually)
- Various fixes to templates
- Updated Terms of Use

## [0.2.2] 09-11-2024

### Changed

- Changed docker image to ubuntu as base

## [0.2.1] 09-11-2024

### Changed

- Pinned `mite_extras` version

## [0.2.0] 09-11-2024

### Changed

- Added submission portal
- Added data parsing/sending

## [0.1.1] 28-10-2024

### Fixed

- Changed docker setup: `prepare_mite_data.py` script must be run beforehand to populate repository
- Separated production (docker) and development (hatch) better

## [0.1.0] 17-10-2024

### Added

- Initialized the repository