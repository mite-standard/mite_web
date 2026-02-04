# Inspired by https://github.com/astral-sh/uv-docker-example/blob/main/standalone.Dockerfile

# First, build the application in the `/app` directory
FROM ghcr.io/astral-sh/uv:bookworm-slim AS builder
ENV UV_COMPILE_BYTECODE=1 UV_LINK_MODE=copy UV_NO_DEV=1

# Configure the Python directory so it is consistent, only use managed version
ENV UV_PYTHON_INSTALL_DIR=/python UV_PYTHON_PREFERENCE=only-managed

# Install Python before the project for caching
RUN uv python install 3.12

# Download and install dependencies
WORKDIR /app
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project

# Install project
COPY . /app
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked

# Download data to be baked in if no data present
ARG DATA=""
ARG EXTRAS=""
RUN set -eux; \
    if [ -d "data" ]; then \
        echo "Using existing data directory - skip download from Zenodo"; \
    elif [ -n "$DATA" ] && [ -n "$EXTRAS" ]; then \
        echo "Downloading data from Zenodo records $DATA and $EXTRAS"; \
        uv run python scripts/prepare_data.py "$DATA" "$EXTRAS" ; \
        uv run python -m scripts.create_db ; \
    else \
        echo "No data directory found and no $DATA and $EXTRAS provided"; \
        echo "   Either populate ./data or pass --build-arg $DATA=..." "$EXTRAS"; \
        exit 1; \
    fi

FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

# Install generic dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    software-properties-common \
    build-essential \
    libxrender1 \
    libxtst6 \
    netcat \
    ncbi-blast+ \
    wget && \
    rm -rf /var/lib/apt/lists/*

# Setup a non-root user
RUN groupadd --system --gid 999 nonroot \
 && useradd --system --gid 999 --uid 999 --create-home nonroot

# Copy the Python version
COPY --from=builder --chown=nonroot:nonroot /python /python

# Copy the application from the builder
COPY --from=builder --chown=nonroot:nonroot /app /app

# Place executables in the environment at the front of the path
ENV PATH="/app/.venv/bin:$PATH"

# Prevent Python from writing .pyc files
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1

# Use the non-root user to run our application
USER nonroot

# Use `/app` as the working directory
WORKDIR /app

# Run the FastAPI application by default
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000", "--proxy-headers", "--forwarded-allow-ips=*"]
