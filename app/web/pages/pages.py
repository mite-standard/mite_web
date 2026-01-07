import os
from pathlib import Path

from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse

from app.config.templates import templates
from app.services.file_handling import load_json
from app.services.items import MiteModel

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


@router.get(
    "/repository/{mite_id}", include_in_schema=False, response_class=HTMLResponse
)
async def repository(mite_id: str, request: Request):
    """Render the static entry pages"""
    try:
        model = MiteModel(mite_id=mite_id)
        return templates.TemplateResponse(
            request=request,
            name="entry.html",
            context={
                "data": load_json(model.html_dir.joinpath(f"{model.mite_id}.json")),
                "next": model.next_id,
                "next_exists": model.next_exists(),
                "prev": model.previous_id,
                "prev_exists": model.previous_exists(),
            },
        )
    except (FileNotFoundError, ValueError):
        return templates.TemplateResponse(
            request=request, name="entry_not_found.html", context={"mite_id": mite_id}
        )

    # TODO: implement testing
