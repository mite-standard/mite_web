mite_web
=========

This repository contains code for web presence of the Minimum Information about a Tailoring Enzyme (MITE) data standard and repository.

The MITE web presence contains functionality to visualize data and receive submissions for new and existing MITE entries.
However, the web presence does not store MITE entries *per se*. The MITE ground truth dataset lives in [mite_data](https://github.com/mite-standard/mite_data).

For more information, see the README of the [MITE-Standard organisation page](https://github.com/mite-standard).

**Nota bene: this application is only tested on Linux and not intended to be run locally. For an online version, see [here](https://mite.bioinformatics.nl/).**

## Installation

### For users

- Install `python3`
- Install `hatch` using one of the methods described [here](https://hatch.pypa.io/1.12/install/)
- Download or clone this repository
- Run `hatch -v env create`

### For developers

- Install `python3`
- Install `hatch` using one of the methods described [here](https://hatch.pypa.io/1.12/install/)
- Download or clone this repository
- Run `hatch -v env create dev`. This will download and install the appropriate Python version and any required packages
- Run `hatch run dev:pre-commit install`. This will set up `pre-commit`

## Quick Start

### For users

- `TBA`

### For developers

- `TBA`
