from http.client import responses

from fastapi import APIRouter, Depends, Request, Form
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
        }
    )
    return response

@router.post("/overview/query", include_in_schema=False, response_class=JSONResponse)
async def overview_query(request: Request):
    response: JSONResponse = JSONResponse(status_code=200, content={"detail": "OK"})
    return response




# @router.post("/overview/query", include_in_schema=False, response_class=HTMLResponse)
# async def overview_query(
#         request: Request, db: Session = Depends(get_db), form: str = Form(...), csrf_protect: CsrfProtect = Depends()
# ):
    # msg = []
    # try:
    #     await csrf_protect.validate_csrf(request)
    #     csrf_token, signed_token = csrf_protect.generate_csrf_tokens()
    #     manager = QueryManager(items={k for k in request.app.state.actives})
    #     # TODO: implement filters
    #
    #     response = templates.TemplateResponse(
    #         request=request,
    #         name="overview.html",
    #         context={
    #             "entries": [v for k, v in request.app.state.actives.items() if k in manager.items],
    #             "headers": request.app.state.table_headers,
    #             "form_vals": request.app.state.form_vals,
    #             "filtered": True,
    #             "form": form
    #         },
    #         csrf_token=csrf_token
    #     )
    #     csrf_protect.set_csrf_cookie(signed_token, response)
    #     return response
    # except Exception as e:
    #     msg.append(
    #         f"An error occurred during search:\n"
    #         f"{e!s}"
    #     )
    #     csrf_token, signed_token = csrf_protect.generate_csrf_tokens()
    #     response = templates.TemplateResponse(
    #         request=request,
    #         name="overview.html",
    #         context={
    #             "entries": [v for v in request.app.state.actives.values()],
    #             "headers": request.app.state.table_headers,
    #             "form_vals": request.app.state.form_vals,
    #             "filtered": False,
    #             "messages": msg
    #         },
    #         csrf_token=csrf_token
    #     )
    #     csrf_protect.set_csrf_cookie(signed_token, response)
    #     return response


# perform query, which keeps track of which entries are removed from set
# once all filters have been performed, the list of actives is returned and used to cut the entries

# TODO: separate route that does the same but streams the result as csv instead.




# @router.post("/overview", include_in_schema=False, response_class=HTMLResponse)
# async def overview(request: Request, db: Session = Depends(get_db), summary: dict = Depends(life)):
