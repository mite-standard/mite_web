from contextlib import asynccontextmanager

from fastapi import FastAPI, Request, staticfiles
from fastapi_csrf_protect import CsrfProtect

from app.api.v1.pages import router_v1
from app.core.config import settings
from app.core.csrf import CsrfSettings
from app.core.templates import configure_templates
from app.db.database import Base, SessionLocal, engine
from app.db.models import Entries
from app.services.file_handling import load_json
from app.web.pages import pages
from app.web.views import overview, views


@asynccontextmanager
async def lifespan(app: FastAPI):
    print("App starting up - DB ready")
    yield
    engine.dispose()
    print("App shutting down - DB disposed")


app = FastAPI(title=settings.app_name, lifespan=lifespan)

app.mount("/data", staticfiles.StaticFiles(directory=settings.data_dir), name="data")
app.mount(
    "/static", staticfiles.StaticFiles(directory=settings.static_dir), name="static"
)
configure_templates(app)


@CsrfProtect.load_config
def get_csrf_config():
    return CsrfSettings()


app.include_router(router_v1)
app.include_router(pages.router)
app.include_router(views.router)
app.include_router(overview.router)


@app.on_event("startup")
def on_startup():
    # TODO: replace with proper seeding
    db = SessionLocal()
    src = settings.data_dir.joinpath("data/MITE0000001.json")
    data = load_json(src)
    entry = Entries(
        accession=data["accession"],
    )
    db.add(entry)
    db.commit()
    db.close()


# TODO: implement on startup database connection, github authentication, S3 connection;
