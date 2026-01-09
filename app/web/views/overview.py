import json

from fastapi import APIRouter, Depends, Form, Request
from fastapi.responses import HTMLResponse, JSONResponse
from sqlalchemy.orm import Session

from app.core.templates import templates
from app.db.crud import QueryManager
from app.db.database import get_db

router = APIRouter(tags=["views"])

# TODO: have a POST route for overview to query the database, consider the streamed csv results
# TODO: complete rework of overview
# TODO: complete rework of overview_retired (needs retired entries dict)


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
        manager = QueryManager(
            items={k for k in request.app.state.actives},
            headers=request.app.state.table_headers,
        )
        forms = dict(await request.form())
        manager.get_entries(forms=forms, db=db)

        # TODO: implement remaining filters

        response = templates.TemplateResponse(
            request=request,
            name="overview.html",
            context={
                "entries": [
                    v
                    for k, v in request.app.state.actives.items()
                    if k in manager.items
                ],
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


@router.post("/overview/download", include_in_schema=False, response_class=HTMLResponse)
async def overview_download(request: Request, db: Session = Depends(get_db)):
    manager = QueryManager(
        items={k for k in request.app.state.actives},
        headers=request.app.state.table_headers,
    )

    form = dict(await request.form())
    forms = json.loads(form["csv_params"])

    manager.get_entries(forms=forms, db=db)

    print([v for k, v in request.app.state.actives.items() if k in manager.items])

    # TODO: continue here by implementing streaming


# perform query, which keeps track of which entries are removed from set
# once all filters have been performed, the list of actives is returned and used to cut the entries

# TODO: separate route that does the same but streams the result as csv instead.
