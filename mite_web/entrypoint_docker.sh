#!/bin/bash

hatch run dev:gunicorn --worker-class gevent --workers 1 "mite_web:create_app()" --bind "0.0.0.0:8004"