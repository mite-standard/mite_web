import json
import logging
import time
import uuid
from collections import defaultdict
from typing import Annotated, Union

from fastapi import APIRouter, Depends, Form, HTTPException, Request
from fastapi.responses import HTMLResponse, StreamingResponse
from github import Github, Repository
from mite_extras import MiteParser
from mite_schema import SchemaManager

from app.auth.basic import get_current_user
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
from app.services.github import (
    add_pr_label,
    create_pr,
    delete_file,
    draft_to_full,
    get_data,
    get_github,
    get_kanban_cached,
    upsert_json_file,
)
from app.services.submission import sign_state, verify_state

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/submission", tags=["views"])


@router.get("/", include_in_schema=False, response_class=HTMLResponse)
async def submission(request: Request, repo: Union[Repository, None] = Depends(get_github)):
    """Get submission page with optional Kanban board"""
    token = sign_state(
        SubmissionState(
            u_id=str(uuid.uuid1()), step="draft", issued=time.time(), role="submitter"
        )
    )
    kanban = None

    if repo:
        kanban = await get_kanban_cached(repo)

    return templates.TemplateResponse(
        request=request,
        name="submission.html",
        context={"kanban": kanban, "token": token},
    )


@router.post("/new", include_in_schema=False, response_class=HTMLResponse)
async def submission_new(
    request: Request,
    form: Annotated[NewDraftForm, Form()],
    repo: Union[Repository, None] = Depends(get_github),
):
    """Create (pull request draft for) new submission"""
    state = verify_state(form.token)
    if state.step != "draft":
        raise HTTPException(400)

    data_model = NewDraftService().parse(form=form)

    if repo:
        await create_pr(repo=repo, branch=state.u_id, data=data_model.data)

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
    repo: Union[Repository, None] = Depends(get_github),
):
    """Create (pull request draft for) existing entry"""
    state = verify_state(form.token)
    if state.step != "draft":
        raise HTTPException(400)

    data_model = ExistDraftService().parse(form=form)

    if repo:
        await create_pr(repo=repo, branch=state.u_id, data=data_model.data)

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
    """Validate form data and redirect to preview based on role"""
    form = await request.form()

    state = verify_state(form["token"])
    if state.step != "preview":
        raise HTTPException(400)

    data = defaultdict(list)
    for key, value in form.multi_items():
        data[key].append(value)
    data = dict(data)

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
                "preview": True,
            },
        )
    else:
        return templates.TemplateResponse(
            request=request,
            name="preview_reviewer.html",
            context={
                "data": model.data.to_html(),
                "data_form": model.data.to_json(),
                "token": token,
                "preview": True,
            },
        )


@router.post("/modified", include_in_schema=False, response_class=HTMLResponse)
async def submission_modified(request: Request):
    """Arrive from preview to do another round of modification"""
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


@router.post("/submit", include_in_schema=False, response_class=HTMLResponse)
async def submission_submit(
    request: Request, repo: Union[Repository, None] = Depends(get_github)
):
    """Submit a submission as a contributor"""
    form = dict(await request.form())

    state = verify_state(form["token"])
    if state.step != "final":
        raise HTTPException(400)

    raw_data = json.loads(form["data_form"])

    if repo:
        await draft_to_full(repo=repo, branch=state.u_id)
        await upsert_json_file(
            repo=repo, branch=state.u_id, data=raw_data, name=state.u_id
        )

    return templates.TemplateResponse(
        request=request,
        name="submission_success.html",
        context={
            "sub_id": state.u_id if repo else None,
        },
    )


@router.post("/download", include_in_schema=False, response_class=StreamingResponse)
async def submission_download(request: Request):
    """Download an entry as JSON"""

    def json_stream(d_in):
        yield json.dumps(d_in, indent=4)

    form = dict(await request.form())
    data = json.loads(form["download-data"])

    return StreamingResponse(
        json_stream(data),
        media_type="application/json",
        headers={"Content-Disposition": 'attachment; filename="mite_entry.json"'},
    )


@router.get("/review/{u_id}", include_in_schema=False, response_class=HTMLResponse)
async def submission_review(
    request: Request,
    u_id: str,
    repo: Union[Repository, None] = Depends(get_github),
    current_user: str = Depends(get_current_user),
):
    """Retrieve data from GitHub API and start review"""
    if not repo:
        HTTPException(400)

    token = sign_state(
        SubmissionState(
            u_id=u_id,
            step="final",
            issued=time.time(),
            role="reviewer",
            reviewer=current_user,
        )
    )

    raw_data = await get_data(repo=repo, branch=u_id)
    model = MiteData(raw_data=raw_data)

    if current_user in raw_data["changelog"][-1]["contributors"]:
        raise HTTPException(
            status_code=400, detail="You can't review your own entries!"
        )

    return templates.TemplateResponse(
        request=request,
        name="preview_reviewer.html",
        context={
            "data": model.data.to_html(),
            "data_form": model.data.to_json(),
            "token": token,
            "preview": True,
        },
    )


@router.post("/review/submit", include_in_schema=False, response_class=HTMLResponse)
async def review_submit(
    request: Request,
    repo: Union[Repository, None] = Depends(get_github),
    current_user: str = Depends(get_current_user),
):
    """Create review and submit to GitHub API"""
    if not repo:
        HTTPException(400)

    form = dict(await request.form())

    state = verify_state(form["token"])
    if state.step != "final" or state.reviewer != current_user:
        raise HTTPException(400)

    raw_data = json.loads(form["data_form"])
    raw_data["changelog"][-1]["reviewers"][0] = state.reviewer
    model = MiteData(raw_data=raw_data)

    if raw_data["accession"] != "MITE9999999":
        await delete_file(repo=repo, branch=state.u_id)

    await upsert_json_file(
        repo=repo,
        branch=state.u_id,
        data=model.data.to_json(),
        name=raw_data["accession"],
    )
    await add_pr_label(repo=repo, branch=state.u_id, label="reviewed")

    return templates.TemplateResponse(
        request=request,
        name="review_success.html",
        context={
            "sub_id": state.u_id if repo else None,
        },
    )
