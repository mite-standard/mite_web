import logging
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from app.api.v1.pages import router_v1
from app.core.config import settings
from app.core.logging import setup_logger
from app.core.shared import load_active, load_form_vals, load_retired, load_table_head
from app.core.templates import configure_templates
from app.db.database import engine
from app.services.github import authenticate_pat
from app.web.pages import download, pages, repository, robots
from app.web.views import overview, pathway, structures, submission

setup_logger()
logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("App starting up - DB ready")

    gh = authenticate_pat()
    app.state.gh = gh

    if gh:
        user = gh.get_user().login
        logger.info(f"Authenticated to GitHub as {user}")
        app.state.repo = gh.get_repo(settings.repo_name)
        logger.info(f"Loaded repo {settings.repo_name}")
    else:
        app.state.repo = None
        logger.error("GitHub authentication failed")

    app.state.actives = load_active()
    app.state.retired = load_retired()
    app.state.table_headers = load_table_head()
    app.state.form_vals = load_form_vals()

    yield

    logger.info("App shutting down.")
    logger.info("DB disposed")
    engine.dispose()
    if app.state.gh:
        app.state.gh.close()
        logger.info("GitHub connection closed.")


app = FastAPI(title=settings.app_name, lifespan=lifespan)

app.mount("/img", StaticFiles(directory=settings.img_dir), name="img")
app.mount("/static", StaticFiles(directory=settings.static_dir), name="static")
configure_templates(app)


app.include_router(router_v1)
app.include_router(pages.router)
app.include_router(overview.router)
app.include_router(pathway.router)
app.include_router(structures.router)
app.include_router(download.router)
app.include_router(robots.router)
app.include_router(repository.router)
app.include_router(submission.router)


# TODO: implement, github authentication
