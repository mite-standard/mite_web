mite_web
=========

# NOTA BENE: FULL REFACTOR IS UNDERWAY. README IS LIKELY OUTDATED


[![DOI](https://zenodo.org/badge/874302233.svg)](https://doi.org/10.5281/zenodo.14933931)

Contents
-----------------
- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Quick Start](#quick-start)
- [Attribution](#attribution)
- [For Developers](#for-developers)

## Overview

**MITE** (Minimum Information about a Tailoring Enzyme) is a community-driven database for the characterization of tailoring enzymes.
These enzymes play crucial roles in the biosynthesis of secondary or specialized metabolites, naturally occurring molecules with strong biological activities, such as antibiotic properties.

This repository manages the code for the [MITE Webpage](https://mite.bioinformatics.nl/).

For more information, visit the [MITE Data Standard Organization page](https://github.com/mite-standard) or read our [publication]( https://doi.org/10.1093/nar/gkaf969).

## Documentation

This repository contains code for the MITE Webpage, allowing to visualize data and receive submissions for new and existing MITE entries.

Note that the web presence **does not store MITE entries** *per se*. The MITE ground truth dataset lives in [mite_data](https://github.com/mite-standard/mite_data).


## System Requirements

**Nota bene: this application is only tested on Linux. For an online version, see [here](https://mite.bioinformatics.nl/).**

### OS Requirements

Local installation was tested on:

- Ubuntu Linux 20.04 and 22.04 (command line)

#### Python dependencies

Dependencies including exact versions are specified in the [pyproject.toml](./pyproject.toml) file.

## Installation Guide

### Deploy Docker locally

- Download or clone this [repository](https://github.com/mite-standard/mite_web)
- Add the `.env` file as indicated below (with at least the **mandatory** params)
- Build the docker image `docker-compose -f docker-compose.yml build --no-cache` (potentially with `sudo`)
- Start the docker `docker-compose -f docker-compose.yml up -d` (potentially with `sudo`)
- Open the application in any browser with the URL http://127.0.0.1:1340/
- To stop the application, run `docker-compose stop` (potentially with `sudo`)
- To take down the database, run `docker-compose down -v --rmi all`

## Attribution

### License

`mite_schema` is an open source tool licensed under the MIT license (see [LICENSE](LICENSE)).

### Publications

See [CITATION.cff](CITATION.cff) or [MITE online](https://mite.bioinformatics.nl/) for information on citing MITE.

### Acknowledgements

This work was supported by the Netherlands Organization for Scientific Research (NWO) KIC grant KICH1.LWV04.21.013.

## For Developers

**Nota bene: since version `1.5.0`, development is only possible using the docker-container**

**Nota bene: regularly update `uv.lock` using `uv lock --upgrade`**

### Development build

- Download or clone this [repository](https://github.com/mite-standard/mite_web)
- Add the `.env` file as indicated below (with at least the **mandatory** params)
- Populate the app by changing into `mite_web` and running `uv run python mite_web/prepare_mite_data.py` (if applicable, remove `.venv`, `data` and `static/img` dirs)
- Build the docker image using `docker-compose build`. This will mount the `mite_web` directory for more convenient file editing (no need to rebuild every time).
- Start the docker image with `docker-compose up -d`. The image will be available at http://127.0.0.1:1340/.
- Changes within the `mite_web` folder will be mirrored inside the docker image but require stopping and restarting the docker container with `docker-compose stop && docker-compose up` (`gunicorn` does not support reloading on change with the current build)
- Changes in the PostgreSQL DB will only be applied if old tables are dropped with `docker-compose down -v`

### Production build

#### First startup

##### 1. Clone the repository

```commandline
git clone https://github.com/mite-standard/mite_web
cd mite_web
```

##### 2. Create configuration file

```commandline
# .env
# mandatory
POSTGRES_PASSWORD=<yoursecurepassword>
POSTGRES_DB=mite_database
# optional
GITHUB_TOKEN=<personal-access-token-classic(scopes: 'admin:public_key', 'gist', 'read:org', 'repo')>
GITHUB_NAME=<gh-acc name>
GITHUB_MAIL=<gh-acc mail>
SECRET_KEY=<your_secret_key>
ONLINE=True
```

##### 3. Build Docker images

```commandline
docker-compose -f docker-compose.yml build --no-cache
```

This builds both `mite_web` and `nginx` images. The mite_web directory is not mounted in production; all code is baked into the image.

##### 4. Start the services

```commandline
docker-compose -f docker-compose.yml up -d
```

##### (5. Remove all services)

```commandline
docker-compose down -v --rmi all
```

#### Update/redeploy

*Nota bene*: make sure to announce the downtime in the [Web App Status](https://github.com/orgs/mite-standard/discussions/5) thread.

##### 1. Save existing dumps (optional)
Preserves links to any open `mite_data` PRs on GitHub

```commandline
docker cp mite_web-mite_web-1:/mite_web/mite_web/dumps .
```

##### 2. Toggle maintenance mode and stop the `mite_web` application 

```commandline
docker exec mite_web-nginx-1 touch /etc/nginx/maintenance.flag
docker exec mite_web-nginx-1 nginx -s reload
docker-compose stop mite_web postgres && docker-compose rm -v postgres
```

This assumes that only the `mite_web` container needs to be updated. `nginx` will continue to run and automatically serve a `maintenance.html` page.

##### 3. Pull the newest release

```commandline
git pull
```

##### 4. Build the `mite_web` image

```commandline
docker-compose -f docker-compose.yml build --no-cache mite_web
```

##### 5. Restart the `mite_web` image and switch off maintenance mode

```commandline
docker-compose -f docker-compose.yml up -d --build --force-recreate postgres && docker-compose -f docker-compose.yml up -d --build --force-recreate mite_web
docker exec mite_web-nginx-1 rm /etc/nginx/maintenance.flag
docker exec mite_web-nginx-1 nginx -s reload
```

##### 6. Transfer existing dumps

```commandline
docker cp ./dumps mite_web-mite_web-1:/mite_web/mite_web && docker cp ./dumps/. mite_web-mite_web-1:/mite_web/mite_web/open_prs
```

##### 7. Perform a test submission

Check if submission system works by performing a test submission. If not, there are possibly issues with the connection to GitHub.