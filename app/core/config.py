from pathlib import Path

from pydantic import Field, model_validator
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """Set up centralized configuration

    Attributes:
        app_name: the name of the app
        static_dir: location of static directory
        data_dir: location of data directory
        img_dir: location of protein images
        github_token: personal-access-token-classic(scopes: 'admin:public_key', 'gist', 'read:org', 'repo'
        github_name: account name
        github_mail: email associated to account
    """

    app_name: str = "Mite Web"
    static_dir: Path = Field(default_factory=lambda: Path("/app/app/static"))
    data_dir: Path = Field(default_factory=lambda: Path("/app/data"))
    img_dir: Path = Field(default_factory=lambda: Path("/app/data/img"))
    github_token: str | None = None
    github_name: str | None = None
    github_mail: str | None = None


settings = Settings()
