from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse

from app.core.templates import templates

router = APIRouter(tags=["pages"])


@router.get("/", include_in_schema=False, response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse(request=request, name="index.html")


@router.get("/termsofuse", include_in_schema=False, response_class=HTMLResponse)
async def termsofuse(request: Request):
    return templates.TemplateResponse(request=request, name="termsofuse.html")


@router.get("/about", include_in_schema=False, response_class=HTMLResponse)
async def about(request: Request):
    return templates.TemplateResponse(request=request, name="about.html")


@router.get("/contact", include_in_schema=False, response_class=HTMLResponse)
async def contact(request: Request):
    return templates.TemplateResponse(request=request, name="contact.html")
