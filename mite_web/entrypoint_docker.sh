#!/bin/bash
echo "Setting up GitHub connection..."
echo "$GITHUB_TOKEN" | gh auth login --with-token --git-protocol https
git config --global credential.helper '!gh auth git-credential'
gh repo clone https://github.com/mite-standard/mite_data.git /mite_web/mite_web/mite_data
git config --global user.name "$GITHUB_NAME"
git config --global user.email "$GITHUB_MAIL"
echo "Completed GitHub connection."

echo "Waiting for postgres..."
while ! nc -z postgres 5432; do
  sleep 0.1
done
echo "PostgreSQL started."

hatch run gunicorn --worker-class gevent --workers 1 "mite_web:create_app()" --bind "0.0.0.0:8004"