mite_web
=========

[![DOI](https://zenodo.org/badge/874302233.svg)](https://doi.org/10.5281/zenodo.14933931)

This repository contains code for web presence of the Minimum Information about a Tailoring Enzyme (MITE) data standard and repository.

The MITE web presence contains functionality to visualize data and receive submissions for new and existing MITE entries.
However, the web presence does not store MITE entries *per se*. The MITE ground truth dataset lives in [mite_data](https://github.com/mite-standard/mite_data).

For more information, see the README of the [MITE-Standard organisation page](https://github.com/mite-standard).

**Nota bene: this application is only tested on Linux and not intended to be run locally. For an online version, see [here](https://mite.bioinformatics.nl/).**

## Deploy Docker locally

- Download or clone this [repository](https://github.com/mite-standard/mite_web)
- Build the docker image `docker-compose build --no-cache` (potentially with `sudo`)
- Start the docker `docker-compose up -d` (potentially with `sudo`)
- Open the application in any browser with the URL http://0.0.0.0:8004/
- To stop the application, run `docker-compose stop` (potentially with `sudo`)

## For developers

*Nota bene: since version `1.5.0`, development is only possible using the docker-container*

### Development build

- Download or clone this [repository](https://github.com/mite-standard/mite_web)
- Create a file `mite_web/instance/config.py` with the content indicated below. Set `Online` to `False`
- Add the `.env` file with content indicated below
- Run the data preparation script from inside the `mite_web` folder using `hatch run python mite_web/prepare_mite_data.py`.
- Remove the hatch environment again with `hatch env remove`
- Build the docker image `docker-compose build`. This will mount the `mite_web` dir for more convenient file editing (no need to rebuild every time).
- Start the docker image with `docker-compose up`
- Changes within the `mite_web` folder will be mirrored inside the docker image but require stopping and restarting the docker container (`gunicorn` does not support reloading on change with the current build)
- Changes in the PostgreSQL DB will only be applied if old tables are dropped with `docker-compose down -v`

### Production build

#### First startup

- Download or clone this [repository](https://github.com/mite-standard/mite_web)
- Create a file `mite_web/instance/config.py` with the content indicated below
- Add the `.env` file with content indicated below
- Build the docker image `docker-compose -f docker-compose.yml build --no-cache` (potentially with `sudo`). Will not mount the `mite_web` dir.
- Start the docker `docker-compose -f docker-compose.yml up -d` (potentially with `sudo`)

#### Update

- Save `dumps` dir (preserves open `mite_data` PR previews) `docker cp mite_web-mite_web-1:/mite_web/mite_web/dumps .`
- To stop the application, run `docker-compose stop` (potentially with `sudo`)
- Take the database down with `docker-compose down -v`
- Pull the newest release
- Build the docker image `docker-compose -f docker-compose.yml build --no-cache` (potentially with `sudo`). Will not mount the `mite_web` dir.
- Start the docker `docker-compose -f docker-compose.yml up -d` (potentially with `sudo`)
- Transfer the `dumps` folder: `docker cp ./dumps mite_web-mite_web-1:/mite_web/mite_web`


#### Config files

`config.py`
```python
SECRET_KEY: str = "your_secret_key"
ONLINE: bool = True
```

`.env`
```commandline
GITHUB_TOKEN=<personal-access-token-classic(scopes: 'admin:public_key', 'gist', 'read:org', 'repo')>
GITHUB_NAME=<gh-acc name>
GITHUB_MAIL=<gh-acc mail>
POSTGRES_PASSWORD=<yoursecurepassword>
POSTGRES_DB=mite_database
```