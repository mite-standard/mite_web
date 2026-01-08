from fastapi import APIRouter, HTTPException, Path, status

from app.services.file_handling import load_json
from app.services.items import MiteListRequest, MiteModel

router_v1 = APIRouter(prefix="/api/v1", tags=["v1"])


@router_v1.post(
    "/mite/",
    summary="Get multiple MITE entries",
    description="Get multiple MITE entry JSONs",
)
async def get_mite_entries(req: MiteListRequest):
    """Return multiple MITE entries"""
    data = []
    for mite_id in req.ids:
        try:
            model = MiteModel(mite_id=mite_id)
            data.append(load_json(model.data_dir.joinpath(f"{model.mite_id}.json")))
        except FileNotFoundError:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail={
                    "error": "File not found",
                    "message": "MITE accession number could not be found.",
                    "requested_path": f"{mite_id}",
                },
            ) from FileNotFoundError
    return data


@router_v1.get(
    "/mite/{mite_acc}",
    summary="Get MITE entry",
    description="Get MITE entry in JSON format",
)
async def get_mite_entry(
    mite_acc: str = Path(
        ...,
        pattern=r"^MITE\d{7}$",
        example="MITE0000001",
    ),
):
    """Return MITE entry as JSON"""
    try:
        model = MiteModel(mite_id=mite_acc)
        return load_json(model.data_dir.joinpath(f"{model.mite_id}.json"))
    except FileNotFoundError:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail={
                "error": "File not found",
                "message": "MITE accession number could not be found.",
                "requested_path": f"/api/v1/mite/{mite_acc}",
            },
        ) from FileNotFoundError
