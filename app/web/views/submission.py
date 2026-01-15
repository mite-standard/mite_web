import logging

from fastapi import APIRouter, Depends, Request
from fastapi.responses import HTMLResponse
from github import Github

from app.core.templates import templates
from app.services.github import get_github, get_kanban_cached

logger = logging.getLogger(__name__)

router = APIRouter(tags=["views"])


@router.get("/submission", include_in_schema=False, response_class=HTMLResponse)
async def submission(request: Request, gh: Github | None = Depends(get_github)):
    if not gh:
        return templates.TemplateResponse(request=request, name="submission.html")

    kanban = await get_kanban_cached(gh)
    return templates.TemplateResponse(
        request=request, name="submission.html", context={"kanban": kanban}
    )
