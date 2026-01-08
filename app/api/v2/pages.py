from fastapi import APIRouter, HTTPException, Path, status

from app.services.file_handling import load_json
from app.services.items import MiteModel

router_v2 = APIRouter(prefix="/api/v2", tags=["v2"])


@router_v2.get(
    "/{mite_acc}",
    summary="MITE Accession ID",
    description="A MITE accession in the format: 'MITE' followed by seven digits (e.g., 'MITE0000001')",
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
                "requested_path": f"/api/v2/{mite_acc}",
            },
        ) from FileNotFoundError
