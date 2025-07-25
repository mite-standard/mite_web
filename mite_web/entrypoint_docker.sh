#!/bin/bash
echo "$GITHUB_TOKEN" | gh auth login --with-token
gh repo clone https://github.com/mite-standard/mite_data.git /mite_web/mite_web/mite_data
git config --global user.name "$GITHUB_NAME"
git config --global user.email "$GITHUB_MAIL"

hatch run gunicorn --worker-class gevent --workers 1 "mite_web:create_app()" --bind "0.0.0.0:8004"