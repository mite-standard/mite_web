from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse

from app.config.templates import templates

router = APIRouter(tags=["views"])


@router.get("/overview", include_in_schema=False, response_class=HTMLResponse)
async def overview(request: Request):
    return templates.TemplateResponse(request=request, name="overview.html")


@router.get("/overview_retired", include_in_schema=False, response_class=HTMLResponse)
async def overview_retired(request: Request):
    return templates.TemplateResponse(request=request, name="overview_retired.html")


@router.get("/submission", include_in_schema=False, response_class=HTMLResponse)
async def submission(request: Request):
    return templates.TemplateResponse(request=request, name="submission.html")


@router.get("/downloads", include_in_schema=False, response_class=HTMLResponse)
async def downloads(request: Request):
    return templates.TemplateResponse(request=request, name="downloads.html")


@router.get("/pathway", include_in_schema=False, response_class=HTMLResponse)
async def pathway(request: Request):
    return templates.TemplateResponse(request=request, name="pathway.html")


@router.get("/peptidesmiles", include_in_schema=False, response_class=HTMLResponse)
async def peptidesmiles(request: Request):
    return templates.TemplateResponse(
        request=request, name="peptidesmiles.html", context={"data": {}}
    )


@router.get("/canonsmiles", include_in_schema=False, response_class=HTMLResponse)
async def canonsmiles(request: Request):
    return templates.TemplateResponse(
        request=request, name="canonsmiles.html", context={"data": {}}
    )
