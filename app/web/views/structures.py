from typing import Annotated

from fastapi import APIRouter, Form, Request
from fastapi.responses import HTMLResponse

from app.core.templates import templates
from app.schemas.structures import GetStructure, PeptideParams, SmilesParams

router = APIRouter(prefix="/structure", tags=["views"])


@router.get("/peptide", include_in_schema=False, response_class=HTMLResponse)
async def peptide(request: Request):
    return templates.TemplateResponse(request=request, name="peptide.html")


@router.post("/peptide/query", include_in_schema=False, response_class=HTMLResponse)
async def peptide_query(request: Request, params: Annotated[PeptideParams, Form()]):
    msg = []
    try:
        runner = GetStructure()
        runner.get_peptide(p=params)

        return templates.TemplateResponse(
            request=request,
            name="peptide.html",
            context={"sequence": params.sequence, "smiles": runner.smiles},
        )
    except Exception as e:
        msg.append(f"An error occurred during peptide smiles generation: {e!s}")
        return templates.TemplateResponse(
            request=request, name="peptide.html", context={"messages": msg}
        )


@router.get("/canonicalize", include_in_schema=False, response_class=HTMLResponse)
async def canonicalize(request: Request):
    return templates.TemplateResponse(request=request, name="canonicalize.html")


@router.post(
    "/canonicalize/query", include_in_schema=False, response_class=HTMLResponse
)
async def canonicalize_query(request: Request, params: Annotated[SmilesParams, Form()]):
    msg = []  # TODO: cleanup message
    try:
        runner = GetStructure()
        runner.get_canonical(p=params)

        return templates.TemplateResponse(
            request=request,
            name="canonicalize.html",
            context={"smiles_in": params.smiles, "smiles_out": runner.smiles},
        )
    except Exception as e:
        msg.append(f"An error occurred during peptide smiles generation: {e!s}")
        return templates.TemplateResponse(
            request=request, name="canonicalize.html", context={"messages": msg}
        )
