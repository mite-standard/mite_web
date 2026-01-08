from typing import Literal

from pydantic import Field, model_validator
from pydantic_settings import BaseSettings


class CsrfSettings(BaseSettings):
    """CSRF configuration for fastapi-csrf-protect

    Attributes:
        secret_key: signs CSRF tokens
        cookie_samesite: defines if cookies can be sent cross-site (not allowed in MITE)
    """

    secret_key: str = Field(..., description="Secret key to sign CSRF tokens")
    environment: Literal["development", "production"] = "development"
    cookie_samesite: Literal["lax", "strict"] = "strict"

    @model_validator(mode="after")
    def validate_security(self):
        if self.environment == "production":
            if self.secret_key == "dev":
                raise RuntimeError("CSRF secret_key must be set in production")
        return self

    @property
    def cookie_secure(self):
        """Cookies sent over https only (production-only)"""
        return self.environment == "production"

    @property
    def cookie_http_only(self):
        """Must be false for CSRF tokens to work"""
        return False

    def dict(self, *args, **kwargs):
        """Accesses properties to materialize them"""
        data = super().dict(*args, **kwargs)
        data["cookie_secure"] = self.cookie_secure
        data["cookie_http_only"] = self.cookie_http_only
        return data
