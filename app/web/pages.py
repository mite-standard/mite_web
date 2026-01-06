from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse

from app.config.templates import templates

router = APIRouter(tags=["pages"])


@router.get("/", include_in_schema=False, response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse(request=request, name="index.html")


@router.get("/tutorial", include_in_schema=False, response_class=HTMLResponse)
async def tutorial(request: Request):
    return templates.TemplateResponse(request=request, name="tutorial.html")


@router.get("/troubleshooting", include_in_schema=False, response_class=HTMLResponse)
async def troubleshooting(request: Request):
    return templates.TemplateResponse(request=request, name="troubleshooting.html")


@router.get("/faqs", include_in_schema=False, response_class=HTMLResponse)
async def faqs(request: Request):
    return templates.TemplateResponse(request=request, name="faqs.html")


@router.get("/termsofuse", include_in_schema=False, response_class=HTMLResponse)
async def termsofuse(request: Request):
    return templates.TemplateResponse(request=request, name="termsofuse.html")


@router.get("/about", include_in_schema=False, response_class=HTMLResponse)
async def about(request: Request):
    return templates.TemplateResponse(request=request, name="about.html")


@router.get("/contact", include_in_schema=False, response_class=HTMLResponse)
async def contact(request: Request):
    return templates.TemplateResponse(request=request, name="contact.html")
