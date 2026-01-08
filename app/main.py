from fastapi import FastAPI, Request, staticfiles
from fastapi_csrf_protect import CsrfProtect

from app.api.v1.pages import router_v1
from app.config.config import CsrfSettings, settings
from app.config.templates import configure_templates
from app.web.pages import pages
from app.web.views import views

app = FastAPI(title=settings.app_name)

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


# TODO: implement on startup database connection, github authentication, S3 connection;
