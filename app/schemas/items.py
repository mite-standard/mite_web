import re
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, computed_field, field_validator, model_validator

MITE_RE = re.compile(r"^MITE(\d{7})$")


class MiteListRequest(BaseModel):
    """Models a list of MITE entries"""

    ids: list[str]


class MiteModel(BaseModel):
    """Model to interact with MITE entries of current release

    Attributes:
        mite_id: mite accession number (=id)
        data_dir: path to data json files
        data_dir: path to html json files
    """

    mite_id: str
    data_dir: ClassVar[Path] = Path("/app/data/data")
    html_dir: ClassVar[Path] = Path("/app/data/html")

    @field_validator("mite_id")
    @classmethod
    def validate_id(cls, v: str) -> str:
        if not MITE_RE.fullmatch(v):
            raise ValueError("Invalid MITE ID format")
        return v

    @model_validator(mode="after")
    def check_exists(self):
        if not self.data_dir.joinpath(f"{self.mite_id}.json").exists():
            raise FileNotFoundError(f"MITE accession does not exist.")
        return self

    @property
    def mite_number(self) -> int:
        """Retrieve the MITE accession number"""
        return int(self.mite_id[4:])

    @computed_field
    @property
    def previous_id(self) -> str | None:
        if self.mite_number <= 1:
            return None
        return f"MITE{self.mite_number - 1:07d}"

    @computed_field
    @property
    def next_id(self) -> str:
        return f"MITE{self.mite_number + 1:07d}"

    def previous_exists(self) -> bool:
        return self.data_dir.joinpath(f"{self.previous_id}.json").exists()

    def next_exists(self) -> bool:
        return self.data_dir.joinpath(f"{self.next_id}.json").exists()
