import json

from fastapi import APIRouter, Depends, Form, Request
from fastapi.responses import HTMLResponse, StreamingResponse
from sqlalchemy.orm import Session

from app.core.templates import templates
from app.db.database import get_db
from app.services.filtering import FilterManager

router = APIRouter(tags=["views"])

# TODO: complete rework of overview


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
async def overview(request: Request):
    response = templates.TemplateResponse(
        "overview.html",
        {
            "request": request,
            "entries": [v for v in request.app.state.actives.values()],
            "headers": request.app.state.table_headers,
            "form_vals": request.app.state.form_vals,
            "filtered": False,
        },
    )
    return response


@router.post("/overview/query", include_in_schema=False, response_class=HTMLResponse)
async def overview_query(request: Request, db: Session = Depends(get_db)):
    msg = []
    try:
        manager = FilterManager(
            accessions={k for k in request.app.state.actives},
            entries=request.app.state.actives,
            headers=request.app.state.table_headers,
        )
        forms = dict(await request.form())
        manager.query_db(forms=forms, db=db)
        manager.query_sequence(forms=forms)
        manager.query_structure(forms=forms)

        # TODO: implement remaining filters

        response = templates.TemplateResponse(
            request=request,
            name="overview.html",
            context={
                "entries": manager.get_entries(),
                "headers": manager.headers,
                "form_vals": request.app.state.form_vals,
                "filtered": True,
                "csv_params": forms,
            },
        )
        return response
    except Exception as e:
        msg.append(f"An error occurred during search:\n" f"{e!s}")
        response = templates.TemplateResponse(
            request=request,
            name="overview.html",
            context={
                "entries": [v for v in request.app.state.actives.values()],
                "headers": request.app.state.table_headers,
                "form_vals": request.app.state.form_vals,
                "filtered": False,
                "messages": msg,
            },
        )
        return response


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

    # TODO: implement rest of filters as above

    return StreamingResponse(
        manager.stream_csv(),
        media_type="text/csv",
        headers={"Content-Disposition": "attachment; filename=query_results.csv"},
    )


# TODO: once all filters have been performed, the list of actives is returned and used to cut the entries
