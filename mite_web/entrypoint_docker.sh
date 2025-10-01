#!/bin/bash

echo "=== Starting Docker entrypoint ==="

echo "Setting up GitHub authentication..."
if [ -z "$GITHUB_TOKEN" ]; then
    echo "ERROR: GITHUB_TOKEN is not set!"
    exit 1
fi
echo "$GITHUB_TOKEN" | gh auth login --with-token --git-protocol https
git config --global credential.helper '!gh auth git-credential'
echo "GitHub authentication configured."

# Check if the mite_data directory already exists
MITE_DATA_DIR="/mite_web/mite_web/mite_data"
if [ -d "$MITE_DATA_DIR" ]; then
    echo "Directory $MITE_DATA_DIR already exists."
    echo "Removing existing directory to ensure a fresh clone..."
    rm -rf "$MITE_DATA_DIR"
    echo "Removed $MITE_DATA_DIR."
else
    echo "Directory $MITE_DATA_DIR does not exist. Ready to clone."
fi

echo "Cloning mite_data repository..."
gh repo clone https://github.com/mite-standard/mite_data.git "$MITE_DATA_DIR"
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to clone the repository."
    exit 1
fi
echo "Repository cloned successfully into $MITE_DATA_DIR."

echo "Configuring Git user..."
git config --global user.name "$GITHUB_NAME"
git config --global user.email "$GITHUB_MAIL"
echo "Git user configured as: $GITHUB_NAME <$GITHUB_MAIL>"

CURRENT_BRANCH=$(git -C "$MITE_DATA_DIR" rev-parse --abbrev-ref HEAD)
CURRENT_COMMIT=$(git -C "$MITE_DATA_DIR" rev-parse --short HEAD)
echo "Current branch: $CURRENT_BRANCH"
echo "Current commit: $CURRENT_COMMIT"

echo "=== GitHub setup complete ==="

echo "Waiting for postgres..."
while ! nc -z postgres 5432; do
  sleep 0.1
done
echo "PostgreSQL started."

uv run gunicorn --worker-class gevent --workers 1 "mite_web:create_app()" --bind "0.0.0.0:8004"