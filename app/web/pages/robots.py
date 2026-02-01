from fastapi import APIRouter
from fastapi.responses import PlainTextResponse

router = APIRouter(tags=["robots"])


@router.get("/robots.txt", response_class=PlainTextResponse)
async def robots_txt():
    return (
        "User - agent: *\n"
        "Disallow: /docs\n"
        "Disallow: /redoc\n"
        "Disallow: /openapi.json\n"
    )
