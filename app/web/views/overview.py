from fastapi import APIRouter, Depends, Request
from fastapi.responses import HTMLResponse
from sqlalchemy.orm import Session

from app.core.templates import templates
from app.crud.entries import get_entires
from app.database import get_db

router = APIRouter(tags=["views"])


@router.get("/overview", include_in_schema=False, response_class=HTMLResponse)
async def overview(request: Request, db: Session = Depends(get_db)):
    print({e.accession for e in get_entires(db=db)})

    # TODO: Complete implementation
    return templates.TemplateResponse(request=request, name="overview.html")
