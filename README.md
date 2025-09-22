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
- Open the application in any browser with the URL http://127.0.0.1:1340/
- To stop the application, run `docker-compose stop` (potentially with `sudo`)

## For developers

*Nota bene: since version `1.5.0`, development is only possible using the docker-container*

### Development build

- Download or clone this [repository](https://github.com/mite-standard/mite_web)
- Create a file `mite_web/instance/config.py` with the content indicated below. Set `Online` to `False`
- Add the `.env` file with content indicated below
- Run the data preparation script from inside the `mite_web` folder using `uv sync && uv run python mite_web/prepare_mite_data.py`.
- Build the docker image `docker-compose build`. This will mount the `mite_web` dir for more convenient file editing (no need to rebuild every time).
- Start the docker image with `docker-compose up`
- Changes within the `mite_web` folder will be mirrored inside the docker image but require stopping and restarting the docker container (`gunicorn` does not support reloading on change with the current build)
- Changes in the PostgreSQL DB will only be applied if old tables are dropped with `docker-compose down -v`

### Production build

#### First startup

1. Clone the repository
```commandline
git clone https://github.com/mite-standard/mite_web
cd mite_web
```
2. Create configuration files
```python
# mite_web/instance/config.py
SECRET_KEY: str = "your_secret_key"
ONLINE: bool = True
```
```commandline
# .env
GITHUB_TOKEN=<personal-access-token-classic(scopes: 'admin:public_key', 'gist', 'read:org', 'repo')>
GITHUB_NAME=<gh-acc name>
GITHUB_MAIL=<gh-acc mail>
POSTGRES_PASSWORD=<yoursecurepassword>
POSTGRES_DB=mite_database
```
3. Build Docker images
```commandline
docker-compose -f docker-compose.yml build --no-cache
```
This builds both `mite_web` and `nginx` images. The mite_web directory is not mounted in production; all code is baked into the image.
4. Start the services
```commandline
docker-compose -f docker-compose.yml up -d
```
5. Remove all services
```commandline
docker-compose down -v --rmi all
```

#### Update/redeploy

*Nota bene*: make sure to announce the downtime in the [Web App Status](https://github.com/orgs/mite-standard/discussions/5) thread.

1. Save existing dumps
```commandline
docker cp mite_web-mite_web-1:/mite_web/mite_web/dumps .
```
optional; preserves open `mite_data` PR previews
2. Stop the `mite_web` application 
```commandline
docker-compose stop mite_web postgres && docker-compose rm -v postgres
```
This assumes that only the `mite_web` container is updated. `nginx` will continue to run and automatically serve a `maintenance.html` page.
3. Pull the newest release
```commandline
git pull
```
4. Build the `mite_web` image
```commandline
docker-compose build --no-cache mite_web
```
5. Restart the `mite_web` image
```commandline
docker compose up -d postgres && docker compose up -d mite_web
```
6. Transfer existing dumps
```commandline
docker cp ./dumps mite_web-mite_web-1:/mite_web/mite_web
```