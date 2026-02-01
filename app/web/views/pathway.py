from typing import Annotated

from fastapi import APIRouter, Form, Request
from fastapi.responses import HTMLResponse, StreamingResponse

from app.core.templates import templates
from app.schemas.pathway import PathwayParams, RunPathway

router = APIRouter(tags=["views"])


@router.get("/pathway/", include_in_schema=False, response_class=HTMLResponse)
async def pathway(request: Request):
    return templates.TemplateResponse(
        request=request,
        name="pathway.html",
    )


@router.post("/pathway/query", include_in_schema=False, response_class=HTMLResponse)
async def pathway_query(request: Request, params: Annotated[PathwayParams, Form()]):
    try:
        runner = RunPathway()
        runner.get_pathway(p=params)

        return templates.TemplateResponse(
            request=request,
            name="pathway.html",
            context={
                "svgs": runner.svgs,
                "report": runner.report,
                "csv_params": params.dump_params(),
            },
        )
    except Exception as e:
        return templates.TemplateResponse(
            request=request,
            name="pathway.html",
            context={
                "messages": f"An error occurred during pathway generation:\n{e!s}"
            },
        )


@router.post("/pathway/csv", include_in_schema=False, response_class=StreamingResponse)
async def pathway_csv(request: Request, params: Annotated[PathwayParams, Form()]):
    runner = RunPathway()
    runner.get_pathway(p=params)
    return StreamingResponse(
        runner.stream_csv(),
        media_type="text/csv",
        headers={"Content-Disposition": "attachment; filename=pathway_results.csv"},
    )
