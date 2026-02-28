#!/usr/bin/env bash
set -euo pipefail

IMAGE="ghcr.io/mite-standard/mite_web"

APP_VERSION="${APP_VERSION:?APP_VERSION not set}"

# mite_data Zenodo: https://doi.org/10.5281/zenodo.13294303
MITE_DATA_VERSION="1.23"
MITE_DATA_RECORD=18806474

# mite_web_extras Zenodo: https://doi.org/10.5281/zenodo.17453501
MITE_WE_VERSION="1.23"
MITE_WE_RECORD=18806707

docker build \
  --build-arg DATA="$MITE_DATA_RECORD" \
  --build-arg EXTRAS="$MITE_WE_RECORD" \
  --tag "$IMAGE:$APP_VERSION" \
  --tag "$IMAGE:$APP_VERSION-data-$MITE_DATA_VERSION-extras-$MITE_WE_VERSION" \
  --tag "$IMAGE:latest" \
  .
