from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse

from app.core.templates import templates
from app.schemas.items import MiteModel
from app.services.file_handling import load_json

router = APIRouter(tags=["pages"])


@router.get(
    "/repository/{mite_id}", include_in_schema=False, response_class=HTMLResponse
)
async def repository(mite_id: str, request: Request):
    """Render the static entry pages"""
    try:
        model = MiteModel(mite_id=mite_id)
        return templates.TemplateResponse(
            request=request,
            name="entry.html",
            context={
                "data": load_json(model.html_dir.joinpath(f"{model.mite_id}.json")),
                "next": model.next_id,
                "next_exists": model.next_exists(),
                "prev": model.previous_id,
                "prev_exists": model.previous_exists(),
            },
        )
    except (FileNotFoundError, ValueError):
        return templates.TemplateResponse(
            request=request, name="entry_not_found.html", context={"mite_id": mite_id}
        )
