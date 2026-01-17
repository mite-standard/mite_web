import re
from datetime import date
from pathlib import Path
from typing import ClassVar, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from app.services import validation
from app.services.file_handling import load_json


class SubmissionState(BaseModel):
    """Define submission state

    u_id: A UUID
    step: Step in submission process
    issued: time in seconds to check for timeout
    role: submitter/reviewer role
    """

    u_id: str
    step: Literal["draft", "preview", "final"]
    issued: float
    role: Literal["submitter", "reviewer"]


class NewDraftForm(BaseModel):
    """New draft entry input form parameters"""

    name: str = Field(description="Enzyme name", max_length=10, min_length=2)
    accession: str = Field(description="UniProt/Genbank ID")
    db: Literal["uniprot", "genpept"]
    submitter: str = Field(description="ORCID or anonymous submission")
    token: str = Field(description="Submission token")

    @field_validator("submitter", mode="before")
    @classmethod
    def normalize_submitter(cls, v: str):
        anonymous = "AAAAAAAAAAAAAAAAAAAAAAAA"
        v = v or anonymous
        validation.validate_submitter(v)
        return v

    @model_validator(mode="after")
    def valid_accession(self):
        if self.db == "uniprot":
            if not validation.is_uniprot(self.accession):
                raise ValueError(f"Not found in UniProtKB/UniParc: {self.accession}")
        else:
            if not validation.is_genpept(self.accession):
                raise ValueError(f"Not found in NCBI GenBank: {self.accession}")
        return self


class NewDraftData(BaseModel):
    """New draft data before data input"""

    data: dict


class NewDraftService:
    """Prepare form data for new draft for storage"""

    @staticmethod
    def get_enzyme(form: NewDraftForm) -> dict:
        if form.db == "uniprot":
            uniprot = form.accession
            genpept = ""
        else:
            genpept = form.accession
            uniprot = ""

        return {
            "name": form.name,
            "databaseIds": {
                "uniprot": uniprot,
                "genpept": genpept,
            },
            "references": ["doi:"],
        }

    def parse(self, form: NewDraftForm) -> NewDraftData:
        return NewDraftData(
            data=dict(
                accession="MITE9999999",
                status="pending",
                changelog=[
                    {
                        "version": "1",
                        "date": date.today().strftime("%Y-%m-%d"),
                        "contributors": [form.submitter],
                        "reviewers": ["BBBBBBBBBBBBBBBBBBBBBBBB"],
                        "comment": "New entry.",
                    }
                ],
                enzyme=self.get_enzyme(form),
                reactions=[
                    {
                        "evidence": {"evidenceCode": [], "references": ["doi:"]},
                        "reactions": [{"products": [""]}],
                    }
                ],
            )
        )


class ExistDraftForm(BaseModel):
    """Existing draft entry input form parameters"""

    mite_id: str = Field(description="MITE Accession ID")
    submitter: str = Field(description="ORCID or anonymous submission")
    token: str = Field(description="Submission token")
    data_dir: ClassVar[Path] = Path("/app/data/data")

    @field_validator("mite_id")
    @classmethod
    def val_id(cls, v: str):
        validation.check_mite_id(v)
        return v

    @field_validator("submitter")
    @classmethod
    def val_submitter(cls, v: str):
        if v == "":
            v = "AAAAAAAAAAAAAAAAAAAAAAAA"
        validation.validate_submitter(v)
        return v

    @model_validator(mode="after")
    def val_exists(self):
        validation.check_exists(self.data_dir.joinpath(f"{self.mite_id}.json"))
        return self


class ExistDraftData(BaseModel):
    """Existing draft data before data input"""

    data: dict


class ExistDraftService:
    """Prepare form data for existing draft for storage"""

    @staticmethod
    def parse(form: ExistDraftForm) -> ExistDraftData:
        """Load existing entry and prepare for form population"""

        data = load_json(form.data_dir.joinpath(f"{form.mite_id}.json"))
        data["status"] = "pending"
        data["changelog"].append(
            {
                "version": len(data["changelog"]) + 1,
                "date": date.today().strftime("%Y-%m-%d"),
                "contributors": [form.submitter],
                "reviewers": ["BBBBBBBBBBBBBBBBBBBBBBBB"],
                "comment": "",
            }
        )
        return ExistDraftData(data=data)
