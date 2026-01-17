from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """Set up centralized configuration

    Attributes:
        static_dir: location of static directory
        data_dir: location of data directory
        img_dir: location of protein images
        env: development or production
        secret: secret key to sign HMAC
        gh_token: GitHub personal-access-token-classic(scopes: 'admin:public_key', 'gist', 'read:org', 'repo'
        gh_name: GitHub account name
        gh_mail: GitHub email associated to account
    """

    static_dir: Path = Field(default_factory=lambda: Path("/app/app/static"))
    data_dir: Path = Field(default_factory=lambda: Path("/app/data"))
    img_dir: Path = Field(default_factory=lambda: Path("/app/data/img"))
    env: Literal["development", "production"] = "development"
    secret: str = "dev"
    gh_token: str | None = None
    gh_name: str | None = None
    gh_mail: str | None = None

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


settings = Settings()
