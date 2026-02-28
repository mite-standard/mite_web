from typing import Annotated

from fastapi import APIRouter, Form, Request
from fastapi.responses import HTMLResponse
from mite_extras import MiteParser

from app.core.templates import templates
from app.schemas.structures import (
    DryrunParams,
    GetStructure,
    PeptideParams,
    SmilesParams,
)

router = APIRouter(prefix="/structure", tags=["views"])


@router.get("/peptide", include_in_schema=False, response_class=HTMLResponse)
async def peptide(request: Request):
    return templates.TemplateResponse(request=request, name="peptide.html")


@router.post("/peptide/query", include_in_schema=False, response_class=HTMLResponse)
async def peptide_query(request: Request, params: Annotated[PeptideParams, Form()]):
    try:
        runner = GetStructure()
        runner.get_peptide(p=params)

        return templates.TemplateResponse(
            request=request,
            name="peptide.html",
            context={"sequence": params.sequence, "smiles": runner.smiles},
        )
    except Exception as e:
        return templates.TemplateResponse(
            request=request,
            name="peptide.html",
            context={
                "messages": f"An error occurred during peptide smiles generation: {e!s}"
            },
        )


@router.get("/canonicalize", include_in_schema=False, response_class=HTMLResponse)
async def canonicalize(request: Request):
    return templates.TemplateResponse(request=request, name="canonicalize.html")


@router.post(
    "/canonicalize/query", include_in_schema=False, response_class=HTMLResponse
)
async def canonicalize_query(request: Request, params: Annotated[SmilesParams, Form()]):
    try:
        runner = GetStructure()
        runner.get_canonical(p=params)

        return templates.TemplateResponse(
            request=request,
            name="canonicalize.html",
            context={"smiles_in": params.smiles, "smiles_out": runner.smiles},
        )
    except Exception as e:
        return templates.TemplateResponse(
            request=request,
            name="canonicalize.html",
            context={
                "messages": f"An error occurred during peptide smiles generation: {e!s}"
            },
        )


@router.get("/flatten", include_in_schema=False, response_class=HTMLResponse)
async def flatten(request: Request):
    return templates.TemplateResponse(request=request, name="flatten.html")


@router.post("/flatten/query", include_in_schema=False, response_class=HTMLResponse)
async def flatten_query(request: Request, params: Annotated[SmilesParams, Form()]):
    try:
        runner = GetStructure()
        runner.strip_chirality(p=params)

        return templates.TemplateResponse(
            request=request,
            name="flatten.html",
            context={"smiles_in": params.smiles, "smiles_out": runner.smiles},
        )
    except Exception as e:
        return templates.TemplateResponse(
            request=request,
            name="flatten.html",
            context={
                "messages": f"An error occurred while removing SMILES stereochemistry: {e!s}"
            },
        )


@router.get("/dryrun", include_in_schema=False, response_class=HTMLResponse)
async def dryrun(request: Request):
    return templates.TemplateResponse(request=request, name="dryrun.html")


@router.post("/dryrun/query", include_in_schema=False, response_class=HTMLResponse)
async def dryrun_query(
    request: Request,
    params: Annotated[DryrunParams, Form()],
):
    data = {
        "accession": "MITE0000000",
        "status": "active",
        "changelog": [
            {
                "version": "1",
                "date": "2024-07-30",
                "contributors": ["0000-0001-6534-6609"],
                "reviewers": ["0000-0001-6534-6609"],
                "comment": "Dummy",
            }
        ],
        "enzyme": {
            "name": "AviG2",
            "references": ["doi:10.1016/j.chembiol.2004.08.016"],
            "databaseIds": {"uniprot": "Q93KW1"},
        },
        "reactions": [
            {
                "tailoring": ["Methylation"],
                "reactionSMARTS": params.r_smarts,
                "reactions": [
                    {
                        "substrate": params.substrate,
                        "products": [params.product],
                        "isIntermediate": False,
                    }
                ],
                "evidence": {
                    "evidenceCode": ["Knock-out studies"],
                    "references": ["doi:10.1016/j.chembiol.2004.08.016"],
                },
            }
        ],
    }

    try:
        parser = MiteParser()
        parser.parse_mite_json(data=data)
        return templates.TemplateResponse(
            request=request,
            name="dryrun.html",
            context={"data": parser.to_html(), "raw": params.model_dump()},
        )
    except Exception as e:
        return templates.TemplateResponse(
            request=request,
            name="dryrun.html",
            context={
                "messages": [f"During data validation, an error has occurred: {e!s}"],
                "raw": params.model_dump(),
            },
        )
