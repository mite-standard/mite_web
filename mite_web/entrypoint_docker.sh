#!/bin/bash

python3 ./prepare_mite_data.py
gunicorn --worker-class gevent --workers 1 "mite_web:create_app()" --bind "0.0.0.0:8002" # TODO(MMZ 21.10.24): add correct IP