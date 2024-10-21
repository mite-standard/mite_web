mite_web
=========

This repository contains code for web presence of the Minimum Information about a Tailoring Enzyme (MITE) data standard and repository.

The MITE web presence contains functionality to visualize data and receive submissions for new and existing MITE entries.
However, the web presence does not store MITE entries *per se*. The MITE ground truth dataset lives in [mite_data](https://github.com/mite-standard/mite_data).

For more information, see the README of the [MITE-Standard organisation page](https://github.com/mite-standard).

**Nota bene: this application is only tested on Linux and not intended to be run locally. For an online version, see [here](https://mite.bioinformatics.nl/).**

## Installation

### For developers

- Install `python3`
- Install `hatch` using one of the methods described [here](https://hatch.pypa.io/1.12/install/)
- Download or clone this repository
- Run `hatch -v env create`. This will download and install the appropriate Python version and any required packages
- Run `hatch run pre-commit install`. This will set up `pre-commit`
- Create a `config.py` file with the content indicated below

## Quick Start

### For developers

Run in dev mode:

- Move into the `mite_web` directory and run `hatch run flask --app mite_web run --debug`

Run as Docker application:

## Config files

For production, a config file must be added to the `mite_web` directory (next to the `pyproject.toml` file)

```python
SECRET_KEY: str
MAIL_USERNAME: str
MAIL_PASSWORD: str
MAIL_DEFAULT_SENDER: str
MAIL_SERVER: str
MAIL_PORT: int
MAIL_USE_TLS: bool
MAIL_USE_SSL: bool
```