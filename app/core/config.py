import base64
import json
from pathlib import Path
from typing import Literal

from pydantic import Field, field_validator, model_validator
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """Set up centralized configuration

    Attributes:
        static_dir: location of static directory
        data_dir: location of data directory
        img_dir: location of protein images
        env: development or production
        reviewers: usernames + pw hashes of reviewers
        secret: secret key to sign HMAC
        gh_token: GitHub personal-access-token-classic(scopes: 'admin:public_key', 'gist', 'read:org', 'repo'
        gh_name: GitHub account name
        gh_mail: GitHub email associated to account
    """

    static_dir: Path = Field(default_factory=lambda: Path("/app/app/static"))
    data_dir: Path = Field(default_factory=lambda: Path("/app/data"))
    img_dir: Path = Field(default_factory=lambda: Path("/app/data/img"))
    env: Literal["development", "production"] = "development"
    reviewers: dict | None = None
    secret: str = "dev"
    gh_token: str | None = None
    gh_name: str | None = None
    gh_mail: str | None = None

    @model_validator(mode="before")
    @classmethod
    def prepare_reviewers(cls, values: dict):
        """Converts base64-utf-8-encoded string into reviewer:pw_hash dict"""
        if not values.get("reviewers"):
            return values

        values["reviewers"] = json.loads(
            base64.b64decode(values["reviewers"]).decode("utf-8")
        )
        return values

    @model_validator(mode="after")
    def check_gh_settings(self):
        gh = [self.gh_token, self.gh_name, self.gh_mail]
        if self.env == "production" and not all(gh):
            raise RuntimeError("One of more GitHub env vars are not set.")
        return self

    @model_validator(mode="after")
    def check_secret(self):
        if self.env == "production" and self.secret == "dev":
            raise RuntimeError("Secret key not set in production")
        return self

    @property
    def repo_name(self) -> str:
        return "mite-standard/mite_data"

    @property
    def app_name(self) -> str:
        return "MITE_Web"

    @property
    def max_age(self) -> int:
        return 60 * 60 * 24 * 7  # 7 days

    @property
    def reviewer_gh_tags(self) -> tuple:
        return ("@mmzdouc",)


settings = Settings()
