import logging
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from app.api.v1.pages import router_v1
from app.core.config import settings
from app.core.logging import setup_logger
from app.core.shared import load_active, load_form_vals, load_table_head
from app.core.templates import configure_templates
from app.db.database import engine
from app.web.pages import pages
from app.web.views import overview, views

setup_logger()
logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("App starting up - DB ready")
    app.state.actives = load_active()
    app.state.table_headers = load_table_head()
    app.state.form_vals = load_form_vals()
    yield
    engine.dispose()
    logger.info("App shutting down - DB disposed")


app = FastAPI(title=settings.app_name, lifespan=lifespan)


app.mount("/data", StaticFiles(directory=settings.data_dir), name="data")
app.mount("/static", StaticFiles(directory=settings.static_dir), name="static")
configure_templates(app)


app.include_router(router_v1)
app.include_router(pages.router)
app.include_router(views.router)
app.include_router(overview.router)


# TODO: implement, github authentication
