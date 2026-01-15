import re
from datetime import date
from typing import Literal

from pydantic import BaseModel, Field, field_validator


class NewDraftForm(BaseModel):
    """Get input for draft form"""

    name: str = Field(description="Enzyme name", max_length=10, min_length=2)
    submitter: str = Field(description="ORCID or anonymous submission")

    @field_validator("submitter")
    @classmethod
    def validate_submitter(cls, v: str):
        if v == "":
            v = "AAAAAAAAAAAAAAAAAAAAAAAA"

        pattern = re.compile(r"^\d{4}-\d{4}-\d{4}-\d{3}[\dX]$|^A{24}$")
        if not re.fullmatch(pattern, v):
            raise ValueError("Not a valid ORCID/not left blank.")

        return v

    def dump(self) -> dict:
        """Dumps content in Github/Form-compatible format"""
        return {
            "accession": "MITE9999999",
            "status": "pending",
            "changelog": [
                {
                    "version": "1",
                    "date": date.today().strftime("%Y-%m-%d"),
                    "contributors": [self.submitter],
                    "reviewers": ["BBBBBBBBBBBBBBBBBBBBBBBB"],
                    "comment": "New entry.",
                }
            ],
            "enzyme": {"name": self.name, "references": ["doi:"]},
            "reactions": [
                {
                    "evidence": {"evidenceCode": [], "references": ["doi:"]},
                    "reactions": [{"products": [""]}],
                }
            ],
        }
