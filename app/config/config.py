import re
from pathlib import Path
from re import Pattern
from typing import ClassVar

from pydantic import Field
from pydantic_settings import BaseSettings

MITE_RE = re.compile(r"^MITE(\d{7})$")


class Settings(BaseSettings):
    """Sets up centralized configuration

    Attributes:
        app_name: the name of the app
        github_token: personal-access-token-classic(scopes: 'admin:public_key', 'gist', 'read:org', 'repo'
        github_name: account name
        github_mail: email associated to account
        debug: development (True), production (False)
        static_dir: location of static directory
        data_dir: location of data directory
    """

    app_name: str = "Mite Web"
    github_token: str | None = None
    github_name: str | None = None
    github_mail: str | None = None
    debug: bool = True
    static_dir: Path = Field(default_factory=lambda: Path("/app/app/static"))
    data_dir: Path = Field(default_factory=lambda: Path("/app/data"))


settings = Settings()
