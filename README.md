mite_web
=========

[![DOI](https://zenodo.org/badge/874302233.svg)](https://doi.org/10.5281/zenodo.14933931)

Contents
-----------------
- [Overview](#overview)
- [Documentation](#documentation)
- [Quickstart Guide](#quickstart-guide)
- [Attribution](#attribution)
- [For Developers](#for-developers)

## Overview

**MITE** (Minimum Information about a Tailoring Enzyme) is a community-driven database for the characterization of tailoring enzymes.
These enzymes play crucial roles in the biosynthesis of secondary or specialized metabolites.
These naturally occurring molecules often show strong biological activities, and many drugs (e.g. antibiotics) derive from them.

This repository manages the code for the [MITE Database](https://mite.bioinformatics.nl/) webapp.

For more information, visit the [MITE Data Standard Organization page](https://github.com/mite-standard) or read our [publication]( https://doi.org/10.1093/nar/gkaf969).

## Documentation

This repository contains code for the MITE Webpage, allowing to visualize data and receive submissions for new and updates for existing entries.

Note that the web presence **does not store MITE entries** *per se*. The MITE ground truth dataset lives in [mite_data](https://github.com/mite-standard/mite_data).

## Quickstart guide

While the [MITE Database](https://mite.bioinformatics.nl/) is primarily intended to be used online, it can also be used offline.

### Installation Guide

*Nota bene: while this application should work on any OS, it has only been tested on Ubuntu Linux 20.04 and 22.04.*
*For an online version, see the [MITE Website](https://mite.bioinformatics.nl/).*

Assuming that Docker is installed:

```commandline
docker run -p 8000:8000 ghcr.io/mite-standard/mite_web:latest
```

You can now use the database running on http://127.0.0.1:8000/.

### Limitations

By default, the application runs in `development` mode. 
Data submission and review are restricted in this mode:

- No submission dashboard
- Newly created entries can only be downloaded, but not submitted to GitHub
- Entries cannot be reviewed.

If entries are to be submitted to the repository, please contact the developers.

## Attribution

### License

`mite_web` is an open source tool licensed under the MIT license (see [LICENSE](LICENSE)).

### Publications

See [CITATION.cff](CITATION.cff) or [MITE online](https://mite.bioinformatics.nl/) for information on citing MITE.

### Acknowledgements

This work was supported by the Netherlands Organization for Scientific Research (NWO) KIC grant KICH1.LWV04.21.013.

## For Developers

### Application logic

*Nota bene: while this application should work on any OS, it has only been tested on Ubuntu Linux 20.04 and 22.04.*

Since `mite_web 2.0.0`, the application follows 12-factor-app principles. 
Data from [mite_data](https://github.com/mite-standard/mite_data) is baked into the image, alongside a read-only SQLite DB for relational querying.
All parameters are provided as environment variables. 
The app itself is stateless, with [mite_data's GitHub repository](https://github.com/mite-standard/mite_data) acting as stateful backend (data transfer via the GitHub API).

### Development build

This build simplifies development by hot reloading (recursively watching directories for changes). 
Data from [mite_data](https://github.com/mite-standard/mite_data) is downloaded once and then injected into the image.
Variables are read from an .env file (see [.env.example](.env.example)).

#### Installation

*Nota bene: assumes that `uv` is installed.*

1. Download and install dependencies
```commandline
git clone git@github.com:mite-standard/mite_data.git
uv sync
uv run pre-commit install
```

2. Run tests
```commandline
uv run pytest
```
Tests will also run via `pre-commit` and GitHub Actions on PRs into main

3. Run script to download data from Zenodo.
```commandline
uv run python scripts/prepare_data.py <mite_data Zenodo record ID> <mite_web_extras Zenodo record ID> # positional args: order must be exact
```
This data is considered dummy data for dev purposes only. 
Make sure to remove it before any (manual) production builds.

4. Run script to build database
```commandline
uv run python -m scripts.create_db
```

5. Run docker compose
```commandline
docker compose -f dev-compose.yml build
docker compose -f dev-compose.yml up --watch
```

Sometimes, the hot reloading doesn't work perfectly (e.g. dependency updates).
In these cases, terminate with `ctrl+c` and rebuild.

#### Deployment Checklist

The main software artifact produced by this repo are Docker containers deposited in the GitHub Container Repository.
These containers are created automatically via GitHub Actions on every new Release.

There are two "types" of releases for Mite Web

1. Only code update
2. Update of code and of baked-in data ([mite_data](https://github.com/mite-standard/mite_data) and [mite_web_extras](https://github.com/mite-standard/mite_web_extras))


##### 1. Only code update

- Branch out
- Perform code updates
- Update version in [`pyproject.toml`](pyproject.toml)
- Describe changes in [`CHANGELOG`](CHANGELOG.md)
- Merge to main
- Create a new release

##### 2. Update of code and baked-in data

- As in 1., with addition of:
- In [build_docker.sh](build_docker.sh), update [mite_data's](https://doi.org/10.5281/zenodo.13294303) `MITE_DATA_VERSION` and `MITE_DATA_RECORD` and [mite_web_extras'](https://doi.org/10.5281/zenodo.17453501) `MITE_WE_VERSION` and `MITE_WE_RECORD`


### Production build

*Nota bene: production build should be exclusively deployed from Docker images (e.g. from [GHCR](ghcr.io/mite-standard/mite_web))*

*Nota bene*: make sure to announce the downtime in the [Web App Status](https://github.com/orgs/mite-standard/discussions/5) thread.

Mite Web requires certain parameters to be run in production mode.
These parameters should be provided as environment variables as specified by the respective platform provider.
An example can be found in [.env.example](.env.example).

If parameters are not provided, the application will automatically start in development mode, where any GitHub API interaction is prevented.

While most parameters are [self-explanatory](.env.example), the reviewer orcids, GitHub tags, and passwords need to be provided as a base64-encoded string.
The easiest way to get this is to create a csv file (see [reviewers.example.csv](reviewers.example.csv)) and run `uv run python scripts/hash_pw.py <input.csv>` on it.

#### Bioinformatics.nl

- Download the content of the [bioinformatics.nl](production/bioinformatics.nl) directory (e.g. with `curl --output <filename> <URL>`)
- Specify the desired container version in [compose.yml](production/bioinformatics.nl/compose.yml)
- Add parameters in an `.env` file (see [.env.example](.env.example)
- Run `docker compose -f compose.yml up`