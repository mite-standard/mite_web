from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles

from app.api.v1.pages import router_v1
from app.config.config import settings
from app.config.templates import configure_templates
from app.web.pages import pages
from app.web.views import views

app = FastAPI(title=settings.app_name)

app.mount("/data", StaticFiles(directory=settings.data_dir), name="data")
app.mount("/static", StaticFiles(directory="/app/app/static"), name="static")

configure_templates(app)

app.include_router(router_v1)
app.include_router(pages.router)
app.include_router(views.router)


# TODO: implement on startup database connection, github authentication, S3 connection;
