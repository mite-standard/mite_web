import json
import re
from datetime import date
from pathlib import Path
from typing import Any, ClassVar, Literal

from mite_extras import MiteParser
from mite_schema import SchemaManager
from pydantic import BaseModel, Field, field_validator, model_validator

from app.services import validation
from app.services.file_handling import load_json


class SubmissionState(BaseModel):
    """Define submission state

    u_id: A UUID
    step: Step in submission process
    issued: time in seconds to check for timeout
    reviewer: authenticated reviewer username
    role: submitter/reviewer role
    """

    u_id: str
    step: Literal["draft", "preview", "final"]
    issued: float
    reviewer: str | None = None
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
                "version": str(len(data["changelog"]) + 1),
                "date": date.today().strftime("%Y-%m-%d"),
                "contributors": [form.submitter],
                "reviewers": ["BBBBBBBBBBBBBBBBBBBBBBBB"],
                "comment": "",
            }
        )
        return ExistDraftData(data=data)


class MiteData(BaseModel):
    """Validated data following mite schema"""

    raw_data: dict
    data: MiteParser

    @model_validator(mode="before")
    @classmethod
    def parse_and_validate(cls, values: dict):
        """Derives MiteParser from values dict and attaches it to the model before creation"""
        parser = MiteParser()
        parser.parse_mite_json(data=values["raw_data"])
        schema_manager = SchemaManager()
        schema_manager.validate_mite(instance=parser.to_json())
        values["data"] = parser
        return values

    @model_validator(mode="after")
    def validate_nr_items(self):
        data = self.data.to_json()
        validation.validate_nr_items(
            container=data["enzyme"]["references"],
            msg="Enzyme: at least one reference DOI needed.",
        )

        for idx, r in enumerate(data["reactions"]):
            validation.validate_nr_items(
                container=r["tailoring"],
                msg="At least one of the checkboxes in 'Tailoring Reaction Controlled Vocabulary' must be checked.",
            )
            validation.validate_nr_items(
                container=r["evidence"]["evidenceCode"],
                msg=f"Reaction {idx +1}: At least one of the checkboxes in 'Experimental Evidence Qualifiers' must be checked.",
            )
            validation.validate_nr_items(
                container=r["evidence"]["references"],
                msg=f"Reaction {idx +1}: At least one reference DOI needed.",
            )

        return self


class MiteService:
    """Data convertion and enrichment"""

    def parse(self, data: dict) -> dict:
        """Convert form data to format compatible with MiteData"""

        return {
            "accession": data["accession"][0],
            "status": data["status"][0],
            "comment": data["comment"][0],
            "changelog": self.changelog(data),
            "enzyme": self.enzyme(data),
            "reactions": self.reactions(data),
        }

    @staticmethod
    def changelog(data: dict) -> list:
        changelog = json.loads(data["changelog"][0])
        changelog[-1]["comment"] = data["changelog-comment"][0]
        return changelog

    def enzyme(self, data: dict) -> dict:
        """Parse out enzyme information"""

        enzyme = {
            "name": data.get("enzyme_name", [""])[0],
            "description": data.get("enzyme_description", [""])[0],
            "references": data.get("enzyme_ref[]", []),
            "databaseIds": {
                "uniprot": self.first_or_none(data, "enzyme_uniprot"),
                "genpept": self.first_or_none(data, "enzyme_genpept"),
                "mibig": self.first_or_none(data, "enzyme_mibig"),
                "wikidata": self.first_or_none(data, "enzyme_wikidata"),
            },
            "cofactors": {
                "organic": data.get(f"enzyme-cofactors-organic-check[]", []),
                "inorganic": data.get(f"enzyme-cofactors-inorganic-check[]", []),
            },
        }

        nr_auxenz = self.get_form_instance(
            data=data, pattern=re.compile(r"auxenzyme\[(\d+)\]")
        )
        if nr_auxenz:
            enzyme["auxiliaryEnzymes"] = []
            for index in nr_auxenz:
                enzyme["auxiliaryEnzymes"].append(
                    {
                        "name": data.get(f"auxenzyme[{index}]name", [""])[0],
                        "description": data.get(f"auxenzyme[{index}]description", [""])[
                            0
                        ],
                        "databaseIds": {
                            "uniprot": data.get(f"auxenzyme[{index}]uniprot", [""])[0],
                            "genpept": data.get(f"auxenzyme[{index}]genpept", [""])[0],
                            "wikidata": data.get(f"auxenzyme[{index}]wikidata", [""])[
                                0
                            ],
                        },
                    }
                )

        return enzyme

    def reactions(self, data: dict) -> list:
        """Parse out reaction information"""
        reactions = {}
        for idx in self.get_form_instance(
            data, pattern=re.compile(r"reaction\[(\d+)\]")
        ):
            reactions[idx] = self.get_instance_known_reactions(data, idx)

        data_mite = []
        for rx, krx in reactions.items():
            data_mite.append(
                {
                    "tailoring": data.get(f"reaction[{rx}]tailoring[]", [""]),
                    "description": data.get(f"reaction[{rx}]description", [""])[0],
                    "reactionSMARTS": data.get(f"reaction[{rx}]smarts", [""])[0],
                    "databaseIds": {
                        "rhea": data.get(f"reaction[{rx}]rhea", [""])[0],
                        "ec": data.get(f"reaction[{rx}]ec", [""])[0],
                    },
                    "evidence": {
                        "evidenceCode": data.get(f"reaction[{rx}]evidencecode[]", [""]),
                        "references": data.get(f"reaction[{rx}]ref[]", [""]),
                    },
                    "reactions": self.get_krxs(data, rx=rx, krx=krx),
                }
            )
        return data_mite

    @staticmethod
    def get_krxs(data: dict, rx: int, krx: list[int]):
        """Parse known reaction dict from form data"""
        examples = []
        for k in krx:
            examples.append(
                {
                    "substrate": data.get(
                        f"reaction[{rx}]knownreaction[{k}]substrate", [""]
                    )[0],
                    "products": data.get(
                        f"reaction[{rx}]knownreaction[{k}]products[]", [""]
                    ),
                    "forbidden_products": data.get(
                        f"reaction[{rx}]knownreaction[{k}]forbiddenproducts[]"
                    ),
                    "isIntermediate": data.get(
                        f"reaction[{rx}]knownreaction[{k}]intermediate", [""]
                    )[0],
                    "description": data.get(
                        f"reaction[{rx}]knownreaction[{k}]description", [""]
                    )[0],
                }
            )
        return examples

    @staticmethod
    def get_form_instance(data: dict, pattern: Any) -> list:
        """Find number of reaction patterns to inform parsing"""
        return sorted(
            {int(match.group(1)) for key in data if (match := pattern.search(key))}
        )

    @staticmethod
    def get_instance_known_reactions(data: dict, idx: int) -> list:
        """Find number of patterns to inform parsing"""
        re_known_r = re.compile(rf"reaction\[{idx}\]knownreaction\[(\d+)\]")
        return sorted(
            {int(match.group(1)) for key in data if (match := re_known_r.search(key))}
        )

    @staticmethod
    def first_or_none(data: dict, key: str) -> str | None:
        val = data.get(key, [""])[0]
        return val if val != "" else None
