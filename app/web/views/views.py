from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse

from app.core.templates import templates

router = APIRouter(tags=["views"])


@router.get("/submission", include_in_schema=False, response_class=HTMLResponse)
async def submission(request: Request):
    # TODO: Complete implementation
    return templates.TemplateResponse(request=request, name="submission.html")


@router.get("/downloads", include_in_schema=False, response_class=HTMLResponse)
async def downloads(request: Request):
    # TODO: Complete implementation
    return templates.TemplateResponse(request=request, name="downloads.html")
