from pathlib import Path

from fastapi import FastAPI
from fastapi.templating import Jinja2Templates

from app.services.file_handling import load_json

templates = Jinja2Templates(directory="/app/app/templates")


def configure_templates(app: FastAPI):
    """Configure globals in template"""
    templates.env.globals["openapi_url"] = app.docs_url
    templates.env.globals["redoc_url"] = app.redoc_url

    version = load_json(Path("/app/data/version_mite_data.json"))
    templates.env.globals["version_mite_data"] = version.get(
        "version_mite_data", "NONE SPECIFIED"
    )
