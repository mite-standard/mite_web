import os
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles

from app.config.templates import configure_templates
from app.web import pages

DATA_DIR = Path(os.environ.get("DATA_DIR", "/app/data"))

app = FastAPI()

app.mount("/data", StaticFiles(directory=DATA_DIR), name="data")
app.mount("/static", StaticFiles(directory="/app/app/static"), name="static")

configure_templates(app)

app.include_router(pages.router)


# TODO: implement on startup database connection, github authentication, S3 connection;
