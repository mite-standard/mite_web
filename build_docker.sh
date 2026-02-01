#!/usr/bin/env bash
set -euo pipefail

IMAGE="ghcr.io/mite-standard/mite_web"

DATA_VERSION="1.21"
DATA_RECORD=17998707

EXTRAS_VERSION="0.5.0"
EXTRAS_RECORD=17999487

docker build \
  --build-arg DATA="$DATA_RECORD" \
  --build-arg EXTRAS="$EXTRAS_RECORD" \
  --tag "$IMAGE:data-$DATA_VERSION-extras-$EXTRAS_VERSION" \
  --tag "$IMAGE" \
  .
