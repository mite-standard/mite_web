import hashlib

from fastapi import APIRouter, Request

from app.core.config import settings

router = APIRouter(prefix="/debug", tags=["pages"])


@router.get("/header", include_in_schema=False)
async def debug(request: Request):
    return dict(request.headers)


@router.get("/scheme", include_in_schema=False)
async def debug_scheme(request: Request):
    return {"scheme": request.url.scheme, "base_url": str(request.base_url)}


@router.get("/secret", include_in_schema=False)
def debug_secret():
    return hashlib.sha256(settings.secret.encode()).hexdigest()
