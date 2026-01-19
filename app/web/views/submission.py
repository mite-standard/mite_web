import json
import logging
import time
import uuid
from collections import defaultdict
from http.client import HTTPException
from typing import Annotated

from fastapi import APIRouter, Depends, Form, Request
from fastapi.responses import HTMLResponse, StreamingResponse
from github import Github
from mite_extras import MiteParser
from mite_schema import SchemaManager

from app.core.templates import templates
from app.schemas.submission import (
    ExistDraftForm,
    ExistDraftService,
    MiteData,
    MiteService,
    NewDraftForm,
    NewDraftService,
    SubmissionState,
)
from app.services.github import create_pr, get_github, get_kanban_cached, push_data
from app.services.submission import sign_state, verify_state

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/submission", tags=["views"])


@router.get("/", include_in_schema=False, response_class=HTMLResponse)
async def submission(request: Request, gh: Github | None = Depends(get_github)):
    token = sign_state(
        SubmissionState(
            u_id=str(uuid.uuid1()), step="draft", issued=time.time(), role="submitter"
        )
    )
    kanban = None

    if gh:
        kanban = await get_kanban_cached(gh)

    return templates.TemplateResponse(
        request=request,
        name="submission.html",
        context={"kanban": kanban, "token": token},
    )


@router.post("/new", include_in_schema=False, response_class=HTMLResponse)
async def submission_new(
    request: Request,
    form: Annotated[NewDraftForm, Form()],
    gh: Github | None = Depends(get_github),
):
    state = verify_state(form.token)
    if state.step != "draft":
        raise HTTPException(400)

    data_model = NewDraftService().parse(form=form)

    if gh:
        create_pr(gh=gh, uuid=state.u_id)
        push_data(gh=gh, uuid=state.u_id, name=state.u_id, data=data_model.data)
        # TODO: implement sending to github via API, use UUID for branch name

    state.step = "preview"
    state.issued = time.time()
    token = sign_state(state)

    return templates.TemplateResponse(
        request=request,
        name="submission_form.html",
        context={
            "data": data_model.data,
            "form_vals": request.app.state.form_vals,
            "token": token,
        },
    )


@router.post("/existing", include_in_schema=False, response_class=HTMLResponse)
async def submission_existing(
    request: Request,
    form: Annotated[ExistDraftForm, Form()],
    gh: Github | None = Depends(get_github),
):
    state = verify_state(form.token)
    if state.step != "draft":
        raise HTTPException(400)

    data_model = ExistDraftService().parse(form=form)

    if gh:
        create_pr(gh=gh, uuid=state.u_id)
        push_data(
            gh=gh,
            uuid=state.u_id,
            name=state.data["enzyme"]["name"],
            data=data_model.data,
        )
        # TODO: implement sending to github via API, use UUID for branch name

    state.step = "preview"
    state.issued = time.time()
    token = sign_state(state)

    return templates.TemplateResponse(
        request=request,
        name="submission_form.html",
        context={
            "data": data_model.data,
            "form_vals": request.app.state.form_vals,
            "token": token,
        },
    )


@router.post("/preview", include_in_schema=False, response_class=HTMLResponse)
async def submission_preview(request: Request):
    form = await request.form()
    data = defaultdict(list)
    for key, value in form.multi_items():
        data[key].append(value)
    data = dict(data)

    state = verify_state(form["token"])
    if state.step != "preview":
        raise HTTPException(400)

    raw_data = MiteService().parse(data=data)
    try:
        model = MiteData(raw_data=raw_data)
    except Exception as e:
        return templates.TemplateResponse(
            request=request,
            name="submission_form.html",
            context={
                "data": raw_data,
                "form_vals": request.app.state.form_vals,
                "token": sign_state(state),
                "messages": [f"During data validation, an error has occurred: {e!s}"],
            },
        )

    state.step = "final"
    state.issued = time.time()
    token = sign_state(state)

    if state.role == "submitter":
        return templates.TemplateResponse(
            request=request,
            name="preview_submitter.html",
            context={
                "data": model.data.to_html(),
                "data_form": model.data.to_json(),
                "token": token,
            },
        )
    else:
        # return reviewer preview template that leads to password-protected confirmation route
        raise RuntimeError("TBA")


@router.post("/modified", include_in_schema=False, response_class=HTMLResponse)
async def submission_modified(request: Request):
    form = dict(await request.form())

    state = verify_state(form["token"])
    if state.step != "final":
        raise HTTPException(400)

    raw_data = json.loads(form["data_form"])

    state.step = "preview"
    state.issued = time.time()
    token = sign_state(state)

    return templates.TemplateResponse(
        request=request,
        name="submission_form.html",
        context={
            "data": raw_data,
            "form_vals": request.app.state.form_vals,
            "token": token,
        },
    )


# post route for final submission


@router.post("/download", include_in_schema=False, response_class=StreamingResponse)
async def submission_download(request: Request):
    def json_stream(d_in):
        yield json.dumps(d_in, indent=4)

    form = dict(await request.form())
    data = json.loads(form["download-data"])

    return StreamingResponse(
        json_stream(data),
        media_type="application/json",
        headers={"Content-Disposition": 'attachment; filename="mite_entry.json"'},
    )
