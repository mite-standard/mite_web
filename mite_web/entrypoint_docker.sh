#!/bin/bash

hatch run gunicorn --worker-class gevent --workers 4 "mite_web:create_app()" --bind "0.0.0.0:8004"