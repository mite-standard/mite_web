from fastapi import APIRouter, Request
from fastapi.responses import FileResponse, HTMLResponse

from app.core.config import settings
from app.core.templates import templates

router = APIRouter(prefix="/download", tags=["pages"])


@router.get("/", include_in_schema=False, response_class=HTMLResponse)
async def download(request: Request):
    return templates.TemplateResponse(request=request, name="download.html")


@router.get("/overview", include_in_schema=False, response_class=FileResponse)
async def download_overview(request: Request):
    return FileResponse(
        path=settings.data_dir.joinpath("summary.csv"),
        filename="summary.csv",
        media_type="application/octet-stream",
    )


@router.get("/smarts", include_in_schema=False, response_class=FileResponse)
async def download_smarts(request: Request):
    return FileResponse(
        path=settings.data_dir.joinpath("download/dump_smarts.csv"),
        filename="dump_smarts.csv",
        media_type="application/octet-stream",
    )


@router.get("/smiles", include_in_schema=False, response_class=FileResponse)
async def download_smiles(request: Request):
    return FileResponse(
        path=settings.data_dir.joinpath("download/dump_smiles.csv"),
        filename="dump_smiles.csv",
        media_type="application/octet-stream",
    )


@router.get("/fasta", include_in_schema=False, response_class=FileResponse)
async def download_fasta(request: Request):
    return FileResponse(
        path=settings.data_dir.joinpath("download/mite_concat.fasta"),
        filename="mite_concat.fasta",
        media_type="application/octet-stream",
    )
