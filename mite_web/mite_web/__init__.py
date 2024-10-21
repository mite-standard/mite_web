"""Application factory of mite_web Flask app.

Copyright (c) 2024-present Mitja Maximilian Zdouc, PhD

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import contextlib
import json
import os
from importlib import metadata
from pathlib import Path

from flask import Flask

from mite_web.routes import bp


def create_app(test_config: dict | None = None) -> Flask:
    """Factory function for Flask app, automatically detected by Flask.

    Arguments:
        test_config: mapping of app configuration for testing purposes

    Returns:
        An instance of the Flask object
    """
    app = Flask(__name__, instance_relative_config=True)
    app = configure_app(app, test_config)

    # TODO(MMZ 21.4.24): check if data repository has been created

    create_instance_path(app)
    register_context_processors(app)
    app.register_blueprint(bp)
    return app


def configure_app(app: Flask, test_config: dict | None = None) -> Flask:
    """Configure the Flask app.

    Arguments:
        app: The Flask app instance
        test_config: mapping of app configuration, can be injected for testing purposes
    """
    app.config["SECRET_KEY"] = "dev"

    if test_config is None:
        app.config.from_pyfile(
            Path(__file__).parent.parent.joinpath("instance/config.py"), silent=True
        )
    else:
        app.config.from_mapping(test_config)
    return app


def verify_data() -> None:
    """Verifies that mite_data was downloaded and made available

    Raises:
        RuntimeError: Data directory not found or empty
    """
    dirpath = Path(__file__).parent.joinpath("data/data_html")
    if not dirpath.exists() or not list(dirpath.iterdir()):
        raise RuntimeError


def create_instance_path(app: Flask):
    """Create the instance path for the Flask app if not available (for testing purposes).

    Arguments:
        app: The Flask app instance
    """
    with contextlib.suppress(OSError):
        os.makedirs(app.instance_path)


def register_context_processors(app: Flask):
    """Register context processors to get access to variables across all pages.

    Arguments:
        app: The Flask app instance
    """

    @app.context_processor
    def metadata_mite_web() -> dict:
        with open(Path(__file__).parent.joinpath("data/version.json")) as infile:
            content = json.load(infile)

        return {
            "version_mite_web": metadata.version("mite_web"),
            "version_mite_data": content.get("version_mite_data"),
        }
