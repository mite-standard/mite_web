import base64
import hashlib
import hmac
import logging
import time

from fastapi import HTTPException
from itsdangerous import BadSignature, SignatureExpired, URLSafeTimedSerializer
from mite_schema import SchemaManager

from app.core.config import settings
from app.schemas.submission import SubmissionState

logger = logging.getLogger(__name__)


def get_state_serializer() -> URLSafeTimedSerializer:
    return URLSafeTimedSerializer(
        secret_key=settings.secret,
        salt="submission-state-v1",
    )


def sign_state(state: SubmissionState) -> str:
    """Sign content"""
    serializer = get_state_serializer()
    return serializer.dumps(state.model_dump())


def verify_state(token: str) -> SubmissionState:
    """Decode and validate signed content."""
    serializer = get_state_serializer()

    try:
        data = serializer.loads(token, max_age=settings.max_age)
        return SubmissionState.model_validate(data)
    except SignatureExpired as e:
        s = "Submission token expired"
        logger.warning(s)
        raise HTTPException(400, detail=s) from e
    except BadSignature as e:
        s = "Submission token signature invalid"
        logger.warning(s)
        raise HTTPException(400, detail=s) from e
    except Exception as e:
        s = f"Error during decoding of signed token: {e!s}"
        logger.warning(s)
        raise HTTPException(400, detail=s) from e


def check_schema(data: dict):
    """Check if data meets schema"""
    SchemaManager().validate_mite(data)


def check_own_entry(user: str, data: dict):
    """Check if reviewer tries to review own entry"""
    if user in data["changelog"][-1]["contributors"]:
        raise HTTPException(
            status_code=400, detail="You must not review your own entries!"
        )


def pop_dummy_reviewer(data: dict) -> dict:
    """Remove dummy reviewer"""
    reviewers = set(data["changelog"][-1]["reviewers"])
    reviewers.discard("BBBBBBBBBBBBBBBBBBBBBBBB")
    data["changelog"][-1]["reviewers"] = list(reviewers)
    return data


def add_reviewer(data: dict, reviewer: str) -> dict:
    """Add reviewer to list"""
    reviewers = set(data["changelog"][-1]["reviewers"])
    reviewers.add(reviewer)
    data["changelog"][-1]["reviewers"] = list(reviewers)
    return data


def set_active(data: dict) -> dict:
    """Set status from pending to active"""
    if data["status"] == "pending":
        data["status"] = "active"
    return data
