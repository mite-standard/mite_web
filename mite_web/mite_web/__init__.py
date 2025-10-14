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
import logging
import os
import sys
from importlib import metadata
from pathlib import Path

import coloredlogs
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_wtf.csrf import CSRFProtect
from mite_schema import SchemaManager

from mite_web.api.mite_api import mite_ns
from mite_web.config.extensions import api, db
from mite_web.routes import bp


def create_app() -> Flask:
    """Factory function for Flask app, automatically detected by Flask.

    Returns:
        An instance of the Flask object
    """
    app = Flask(__name__, instance_relative_config=True)
    app = configure_app(app)

    db.init_app(app)
    from mite_web.seed import seed_data

    with app.app_context():
        db.create_all()
        seed_data()

    app.url_map.strict_slashes = False
    verify_data(app)

    register_context_processors(app)
    app.register_blueprint(bp)

    api.init_app(app)
    api.add_namespace(mite_ns, path="/api/v1/mite")

    return app


def configure_app(app: Flask) -> Flask:
    """Configure the Flask app.

    Arguments:
        app: The Flask app instance
    """
    app = config_logger(app)
    app = set_paths(app)

    app.config["SECRET_KEY"] = os.getenv("SECRET_KEY", "dev")
    if app.config["SECRET_KEY"] == "dev":
        app.logger.critical("INSECURE DEV MODE: DO NOT DEPLOY TO PRODUCTION!")

    app.config["ONLINE"] = os.getenv("ONLINE", "False")
    if app.config["ONLINE"] == "False":
        app.logger.critical("OFFLINE MODE: WILL BLOCK PRS TO GITHUB!")

    app.config["DATA_DUMPS"].mkdir(parents=True, exist_ok=True)
    app.config["QUERIES"].mkdir(parents=True, exist_ok=True)
    app.config["OPEN_PRS"].mkdir(parents=True, exist_ok=True)

    app = diff_active_retired(app)
    app = populate_form_data(app)

    app.config["SQLALCHEMY_DATABASE_URI"] = os.getenv("DATABASE_URL")
    app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False

    csrf = CSRFProtect()
    csrf.init_app(app)

    return app


def set_paths(app: Flask) -> Flask:
    """Set all paths required by the program

    Args:
        app: The Flask app

    Returns:
        The Flask app with defined paths
    """
    app.config["DATA_HTML"] = Path(__file__).parent.joinpath("data/data_html")
    app.config["DATA_JSON"] = Path(__file__).parent.joinpath("data/data")
    app.config["DATA_DUMPS"] = Path(__file__).parent.joinpath("dumps")
    app.config["OPEN_PRS"] = Path(__file__).parent.joinpath("open_prs")
    app.config["DOWNLOAD"] = Path(__file__).parent.joinpath("data/download")
    app.config["QUERIES"] = Path(__file__).parent.joinpath("queries")
    app.config["DATA_IMG"] = Path(__file__).parent.joinpath("static/img")
    app.config["DATA_SUMMARY"] = Path(__file__).parent.joinpath("data/summary.json")
    app.config["MITE_DATA"] = Path(__file__).parent.joinpath("mite_data")
    return app


def diff_active_retired(app: Flask) -> Flask:
    """Populate separate lists for active and retired entries to show in overview page

    Arguments:
        app: the Flask app as reference

    Returns:
        The unmodified Flask app
    """
    with open(app.config["DATA_SUMMARY"]) as infile:
        summary = json.load(infile)
        app.config["SUMMARY_ACTIVE"] = {
            key: val
            for key, val in summary["entries"].items()
            if val.get("status_plain") == "active"
        }
        app.config["ACCESSIONS_ACTIVE"] = {acc for acc in app.config["SUMMARY_ACTIVE"]}
        app.config["SUMMARY_RETIRED"] = {
            key: val
            for key, val in summary["entries"].items()
            if val.get("status_plain") == "retired"
        }
    return app


def populate_form_data(app: Flask) -> Flask:
    """Retrieve dropdown lists to populate forms

    Args:
        app: The Flask aap

    Returns:
        The unmodified app
    """
    with open(SchemaManager().entry) as f:
        schema = json.load(f)
        app.config["FORM_VALS"] = {
            "evidence": schema["$defs"]["evidence"]["enum"],
            "tailoring": schema["$defs"]["tailoringFunction"]["enum"],
            "inorganic": schema["$defs"]["inorganic"]["enum"],
            "organic": schema["$defs"]["organic"]["enum"],
        }
    return app


def config_logger(app: Flask) -> Flask:
    """Set up a named logger with nice formatting and attach to app

    Args:
        app: The Flask app

    Returns:
        The Flask app with attached logger
    """
    for handler in app.logger.handlers[:]:
        app.logger.removeHandler(handler)

    logger = logging.getLogger("mite_web")
    logger.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(
        coloredlogs.ColoredFormatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
    )

    file_handler = logging.FileHandler(
        Path(__file__).parent.joinpath("app.log"),
        mode="w",
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )

    app.logger.addHandler(console_handler)
    app.logger.addHandler(file_handler)
    return app


def verify_data(app: Flask):
    """Verifies that mite_data was downloaded and made available

    Arguments:
        app: The Flask app

    Raises:
        RuntimeError: Data not found or empty
    """
    if not app.config["DATA_SUMMARY"].exists():
        message = f"Could not find file '{app.config["DATA_SUMMARY"].resolve()}' - did you run the 'prepare_mite_data.py' script?"
        app.logger.critical(message)
        raise RuntimeError(message)

    if not app.config["DATA_HTML"].exists() or not list(
        app.config["DATA_HTML"].iterdir()
    ):
        message = f"Could not find folder '{app.config["DATA_HTML"].resolve()}' - did you run the 'prepare_mite_data.py' script?"
        app.logger.critical(message)
        raise RuntimeError(message)

    if not app.config["DATA_IMG"].exists() or not list(
        app.config["DATA_IMG"].iterdir()
    ):
        message = f"Could not find folder '{app.config["DATA_IMG"].resolve()}' - did you run the 'prepare_mite_data.py' script?"
        app.logger.critical(message)
        raise RuntimeError(message)


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
