from collections import defaultdict

from fastapi import APIRouter, Depends, Request
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import HTMLResponse
from github import Github

from app.core.config import settings
from app.core.templates import templates
from app.services.github import fake_pull, get_github, process_pulls

router = APIRouter(tags=["views"])


@router.get("/submission", include_in_schema=False, response_class=HTMLResponse)
async def submission(request: Request, gh: Github | None = Depends(get_github)):
    def fetch_pulls():
        repo = gh.get_repo(settings.repo_name)
        pulls = repo.get_pulls(state="open")
        return process_pulls(pulls)

    if not gh:
        # TODO: replace with data = None
        return templates.TemplateResponse(
            request=request, name="submission.html", context={"kanban": fake_pull()}
        )

    kanban = await run_in_threadpool(fetch_pulls)
    return templates.TemplateResponse(
        request=request, name="submission.html", context={"kanban": kanban}
    )
