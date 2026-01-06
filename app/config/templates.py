from fastapi import FastAPI
from fastapi.templating import Jinja2Templates

templates = Jinja2Templates(directory="/app/app/templates")


def configure_templates(app: FastAPI):
    """Configure globals in template"""
    templates.env.globals["openapi_url"] = app.docs_url
    templates.env.globals["redoc_url"] = app.redoc_url
