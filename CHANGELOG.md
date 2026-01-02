# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.7.4] 02-01-2026

### Fixed

- Adjusted number of workers back to 1 due to errors in database seeding

## [1.7.3] 07-12-2025

### Changed

- Moved validations of MIBiG and RHEA IDs to `mite_data`

## [1.7.2] 01-11-2025

### Fixed

- Changed routing from `GET` to `POST` to prevent `GET` induced changes in state
- Fixed failing "delete" on in silico synthesis page

## [1.7.1] 28-10-2025

### Fixed

- Adjusted logic for reviewer-modification of entries


## [1.7.0] 27-10-2025

### Changed

- Simplified docker build - auxiliary files are downloaded
- Reworked submission portal - show pending and submitted entries

## [1.6.6] 15-10-2025

### Changed

- Prevent already existing PRs to be re-created

## [1.6.5] 14-10-2025

### Added

- File lock for GitHub PRs
- Overview of currently "pending" entries on GitHub

### Changed

- Rework offline installation
- 'Download file' logic

## [1.6.4] 30-09-2025

### Changed

- Fixes for pushing submissions to GitHub (old `mite_data` folder persisted)

## [1.6.3] 30-09-2025

### Changed

- Added an explicit `maintenance.html` page to catch errors
- Updated `CITATION.cff` for publication in *NAR*

## [1.6.2] 22-09-2025

### Changed

- Replaced `hatch` with `uv` for more reproducible builds

## [1.6.1] 29-08-2025

### Fixed

- Fixed issue where non-filled query rule would prevent search submission

## [1.6.0] 29-08-2025

### Added

- Added Overview page for retired entries
- Implemented in silico biosynthesis page
- Embedded tutorial videos

### Changed

- Removed retired entries from ovrview table
- Updated website descriptions

### Fixed

- Implemented escaping of problematic characters for input fields

## [1.5.0] 10-08-2025

### Changed

- Reworked search page, added complex search
- Added reaction SMARTS diff fingerprint search
- Implemented search query system

## [1.4.2] 26-07-2025

### Fixed

- Added formatting to changelog strings

## [1.4.2] 26-07-2025

### Fixed

- Updated dependencies
- Removed old "review" route

## [1.4.1] 25-07-2025

### Fixed

- Added `--fix-missing` option to Dockerfile for more robust builds

### Changed

## [1.4.0] 25-07-2025

### Added

- Additional source info in overview table (e.g. Domain, Kingdom, ...)
- Download of overview table as .csv
- Forms: added cofactor selection and wikidata cross-link for enzyme info

### Changed

- Forms: options for tailoring, evidence etc. now directly from `mite_schema` data model
- Removed redundant user input sanitation (already performed by `mite_extras`)
- Reworked data submission - changed to direct pull requests to GitHub
- Removed email functionality
- Updated api to be properly url-versioned

## [1.3.1] 02-06-2025

### Fixed

- Frontend fixes
- Gunicorn vulnerability fix

## [1.3.0] 04-03-2025

### Added

- Search for reaction SMARTS

### Changed

- Reworked search functionality on overview page

## [1.2.5] 02-03-2025

### Changed

- Added "scheduled maintenance" to navbar

## [1.2.4] 28-02-2025

### Changed

- Updated metadata to refer to community-level collaboration documents

## [1.2.3] 27-02-2025

### Added

- Added Reviewer page to the Navbar

## [1.2.2] 26-02-2025

### Added

- Implement a page to visualize a MITE JSON file (for reviewing purposes)

### Changed

- MIBiG URL resolve via Bioregistry

## [1.2.1] 14-02-2025

### Bugfix

- Fixed Rhea validation behavior

## [1.2.0] 13-02-2025

### Added

- Added API documentation under `/api/`

### Changed

- Reworked UI (Navbar)
- Changed templates

## [1.1.2] 02-01-2025

### Added

- Added concatenated protein FASTA files to downloads

### Bugfix

- Fixed broken hyperlinks in substructure/blast search results

### Changed

- DatabaseIds crosslinking: differentiate between UniprotKB and UniParc in link-creation (UniParcID starts with "UPI")
- Replaced generic MITE contributor "AAAAA..." with a human-readable explanation
- Moved enzyme-structure visualization with PyMol from `mite_data` to `mite_web`
- Moved MITE BLAST-DB generation from `mite_data` to `mite_web`

## [1.1.1] 09-12-2024

### Added

- Implemented an API under `/api/<mite_accession>` to serve MITE JSON files

### Changed

- Improved search functionality for substructure

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