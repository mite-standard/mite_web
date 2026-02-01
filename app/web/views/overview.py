import json
from urllib.parse import urlencode

from fastapi import APIRouter, Depends, Request
from fastapi.responses import HTMLResponse, RedirectResponse, StreamingResponse
from sqlalchemy.orm import Session

from app.core.templates import templates
from app.db.database import get_db
from app.services.filtering import FilterManager

router = APIRouter(tags=["views"])


@router.get("/retired", include_in_schema=False, response_class=HTMLResponse)
async def retired(request: Request):
    response = templates.TemplateResponse(
        "overview_retired.html",
        {
            "request": request,
            "entries": [v for v in request.app.state.retired.values()],
            "headers": request.app.state.table_headers,
        },
    )
    return response


@router.get("/overview", include_in_schema=False, response_class=HTMLResponse)
async def overview(request: Request, db: Session = Depends(get_db)):
    params = dict(request.query_params)

    filtered = bool(params)
    entries = [v for v in request.app.state.actives.values()]
    headers = request.app.state.table_headers
    messages = []

    if filtered:
        try:
            manager = FilterManager(
                accessions={k for k in request.app.state.actives},
                entries=request.app.state.actives,
                headers=headers,
            )

            manager.query_db(forms=params, db=db)
            manager.query_sequence(forms=params)
            manager.query_structure(forms=params)
            manager.query_reaction(forms=params)

            entries = manager.get_entries()
            headers = manager.headers
        except Exception as e:
            messages.append(f"An error occurred during search:\n{e!s}")

    return templates.TemplateResponse(
        request=request,
        name="overview.html",
        context={
            "entries": entries,
            "headers": headers,
            "form_vals": request.app.state.form_vals,
            "filtered": filtered,
            "csv_params": params,
            "messages": messages,
        },
    )


@router.post(
    "/overview/query", include_in_schema=False, response_class=RedirectResponse
)
async def overview_query(request: Request):
    forms = dict(await request.form())

    return RedirectResponse(
        url=f"/overview?{urlencode(forms)}",
        status_code=303,
    )


@router.post("/overview/csv", include_in_schema=False, response_class=StreamingResponse)
async def overview_csv(request: Request, db: Session = Depends(get_db)):
    manager = FilterManager(
        accessions={k for k in request.app.state.actives},
        entries=request.app.state.actives,
        headers=request.app.state.table_headers,
    )
    form = dict(await request.form())
    forms = json.loads(form["csv_params"])

    manager.query_db(forms=forms, db=db)
    manager.query_sequence(forms=forms)
    manager.query_structure(forms=forms)
    manager.query_reaction(forms=forms)
    return StreamingResponse(
        manager.stream_csv(),
        media_type="text/csv",
        headers={"Content-Disposition": "attachment; filename=query_results.csv"},
    )
