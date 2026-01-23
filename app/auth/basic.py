import secrets

from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPBasic, HTTPBasicCredentials

from app.auth.passwords import verify_pw
from app.core.config import settings

security = HTTPBasic()


def get_current_user(credentials: HTTPBasicCredentials = Depends(security)) -> str:
    """Authenticate reviewer"""
    exception = HTTPException(
        401,
        detail="Invalid credentials",
        headers={"WWW-Authenticate": "Basic"},
    )

    username = credentials.username
    password = credentials.password

    stored_hash = settings.reviewers.get(username)

    if not stored_hash:
        verify_pw(
            plain=password,
            hashed="$2b$12$rbOOwEBSoM3Nh0DepFm5V.br8QYSRboi4wvhrcfQL5QIqSfS44.Ou",
        )
        raise exception

    if not verify_pw(plain=password, hashed=stored_hash):
        raise exception

    return username
