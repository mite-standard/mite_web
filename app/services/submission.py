import base64
import hashlib
import hmac
import time

from app.core.config import settings
from app.schemas.submission import SubmissionState


def sign_state(state: SubmissionState) -> str:
    """Sign content with HMAC logic"""
    payload = state.model_dump_json().encode()
    sig = hmac.new(
        settings.secret.encode(),
        payload,
        hashlib.sha256,
    ).digest()

    return base64.urlsafe_b64encode(payload + b"." + sig).decode()


def verify_state(token: str) -> SubmissionState:
    """Decode signed content."""
    raw = base64.urlsafe_b64decode(token.encode())
    payload, sig = raw.rsplit(b".", 1)

    expected = hmac.new(
        settings.secret.encode(),
        payload,
        hashlib.sha256,
    ).digest()

    if not hmac.compare_digest(sig, expected):
        raise ValueError("Invalid signature")

    state = SubmissionState.model_validate_json(payload)

    if time.time() - state.issued > settings.max_age:
        raise ValueError("Token expired")

    return state
