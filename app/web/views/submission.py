import logging
import uuid
from typing import Annotated

from fastapi import APIRouter, Depends, Form, Request
from fastapi.responses import HTMLResponse
from github import Github

from app.core.templates import templates
from app.schemas.submission import NewDraftForm
from app.services.github import create_pr, get_github, get_kanban_cached, push_data

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/submission", tags=["views"])


@router.get("/", include_in_schema=False, response_class=HTMLResponse)
async def submission(request: Request, gh: Github | None = Depends(get_github)):
    if not gh:
        return templates.TemplateResponse(request=request, name="submission.html")

    kanban = await get_kanban_cached(gh)
    return templates.TemplateResponse(
        request=request,
        name="submission.html",
        context={"kanban": kanban, "csrftoken": request.state.csrftoken},
    )


@router.post("/new", include_in_schema=False, response_class=HTMLResponse)
async def submission_new(
    request: Request,
    form: Annotated[NewDraftForm, Form()],
    gh: Github | None = Depends(get_github),
):
    try:
        u_id = str(uuid.uuid1())
        if gh:
            create_pr(gh=gh, uuid=u_id)
            push_data(gh=gh, uuid=u_id, name=u_id, data=form.dump())
            # TODO: implement sending to github via API, use UUID for branch name

        return templates.TemplateResponse(
            request=request,
            name="submission_form.html",
            context={
                "data": form.dump(),
                "u_id": u_id,
                "form_vals": request.app.state.form_vals,
            },
        )
    except Exception as e:
        msg = [f"{e!s}"]
        if not gh:
            return templates.TemplateResponse(
                request=request, name="submission.html", context={"messages": msg}
            )

        kanban = await get_kanban_cached(gh)
        return templates.TemplateResponse(
            request=request,
            name="submission.html",
            context={"kanban": kanban, "messages": msg},
        )

    # get input values
    # check if setting given
    # if yes, send to function to deposit to github via API
    # redirect to form page with pre-filled values
    # check for empty dummy field to prevent bos submissions


# TODO: post route for data validation and rending of final page, or return loop to a submission_existing with the correct name
