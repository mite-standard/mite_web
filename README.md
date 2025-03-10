mite_web
=========

[![DOI](https://zenodo.org/badge/874302233.svg)](https://doi.org/10.5281/zenodo.14933931)

This repository contains code for web presence of the Minimum Information about a Tailoring Enzyme (MITE) data standard and repository.

The MITE web presence contains functionality to visualize data and receive submissions for new and existing MITE entries.
However, the web presence does not store MITE entries *per se*. The MITE ground truth dataset lives in [mite_data](https://github.com/mite-standard/mite_data).

For more information, see the README of the [MITE-Standard organisation page](https://github.com/mite-standard).

**Nota bene: this application is only tested on Linux and not intended to be run locally. For an online version, see [here](https://mite.bioinformatics.nl/).**

## Deploy Docker

### Quick start

- Download or clone this [repository](https://github.com/mite-standard/mite_web)
- Create a file `mite_web/instance/config.py` with the content indicated below
- Build the docker image `docker-compose build --no-cache` (potentially with `sudo`)
- Start the docker `docker-compose up -d` (potentially with `sudo`)
- Open the application in any browser with the URL http://0.0.0.0:8004/
- To stop the application, run `docker-compose stop` (potentially with `sudo`)

## For developers

### Install and quick start

*Nota bene: This installation will only work on (Ubuntu) Linux.*

- Install `python3`
- Install `hatch` using one of the methods described [here](https://hatch.pypa.io/1.12/install/)
- Download or clone this [repository](https://github.com/mite-standard/mite_web) and change into `mite_dev` (where the `pyproject.toml` file resides)
- Run `hatch -v env create dev`
- Run `hatch run dev:pre-commit install`. This will set up `pre-commit`
- Install BLAST+ by running `sudo apt-get install ncbi-blast+`
- Install Pymol by running `sudo apt-get install pymol`
- Run the script `hatch run dev:python mite_web/prepare_mite_data.py` to populate the application with data
- Run the script `./run_pymol.sh` to populate the application with data
- Move into the `mite_web` directory and run `hatch run dev:flask --app mite_web run --debug`

## Config-file

```python
SECRET_KEY: str = "your_secret_key"
ONLINE: bool = False
MAIL_TARGET: str = "email_target"
MAIL_DEFAULT_SENDER: str = "sender"
MAIL_SERVER: str = "server"
MAIL_PORT: int = "port"
MAIL_USE_TLS: bool = False
MAIL_USE_SSL: bool = False
```
