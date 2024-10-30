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

import json
from importlib import metadata
from pathlib import Path

import wtforms_json
from flask import Flask

from mite_web.routes import bp


def create_app() -> Flask:
    """Factory function for Flask app, automatically detected by Flask.

    Returns:
        An instance of the Flask object
    """
    app = Flask(__name__, instance_relative_config=True)
    app = configure_app(app)
    verify_data()
    register_context_processors(app)
    app.register_blueprint(bp)
    return app


def configure_app(app: Flask) -> Flask:
    """Configure the Flask app.

    Arguments:
        app: The Flask app instance
    """
    app.config["SECRET_KEY"] = "dev"

    wtforms_json.init()

    config_file = Path(__file__).parent.parent.joinpath("instance/config.py")
    if config_file.exists():
        app.config.from_pyfile(config_file)
        print("Successfully loaded configuration from 'config.py'.")
    else:
        print("WARNING: No 'config.py' file found. Default to dev settings.")
        print("WARNING: INSECURE DEV MODE: DO NOT DEPLOY TO PRODUCTION!")

    return app


def verify_data() -> None:
    """Verifies that mite_data was downloaded and made available

    Raises:
        RuntimeError: Data directory not found or empty
    """
    dirpath = Path(__file__).parent.joinpath("data/data_html")
    if not dirpath.exists() or not list(dirpath.iterdir()):
        raise RuntimeError(
            f"Could not find folder '{dirpath.resolve()}' - did you run the 'prepare_mite_data.py' script?"
        )

    imgpath = Path(__file__).parent.joinpath("static/img")
    if not imgpath.exists() or not list(imgpath.iterdir()):
        raise RuntimeError(
            f"Could not find folder '{imgpath.resolve()}' - did you run the 'prepare_mite_data.py' script?"
        )


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
