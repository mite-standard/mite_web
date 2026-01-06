from fastapi import APIRouter, Request
from fastapi.responses import JSONResponse

router_v1 = APIRouter(prefix="/api/v1/mite", tags=["deprecated"])


@router_v1.get(
    "/{string:path}",
    summary="MITE Accession ID",
    description="A MITE accession in the format: 'MITE' followed by seven digits (e.g., 'MITE0000001')",
)
async def deprecated_v1(mite_acc: str, request: Request):
    return JSONResponse(
        status_code=410,
        content={
            "error": "API v1 is deprecated",
            "message": "Please migrate to /api/v2",
            "requested_path": f"/api/v1/{mite_acc}",
        },
    )
