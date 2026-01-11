import csv
import io
import json
import re
import subprocess
import tempfile
from typing import Any, ClassVar

from Bio import Blast
from pydantic import BaseModel, Field, field_validator
from sqlalchemy import ColumnElement, and_, inspect, or_
from sqlalchemy.orm import Session

from app.core.config import settings
from app.db.models import Entry


class BlastManager(BaseModel):
    """Organizes querying for a protein sequence"""

    sequence: str = Field(..., description="Protein sequence")
    e_val: int = Field(10, gt=0, description="E-value for BLASTp search")

    @field_validator("sequence", mode="before")
    @classmethod
    def clean_sequence(cls, v: str) -> str:
        if not isinstance(v, str):
            raise TypeError("Sequence must be a string")

        v = v.replace("\n", "").upper()
        v = re.sub(r"\s+", "", v, flags=re.UNICODE)

        if not v:
            raise ValueError("Sequence is empty")

        if not re.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", v):
            raise ValueError("Sequence contains invalid amino acid characters")

        return v

    def run_blastp(self) -> dict:
        """Run BLAST against MITE BLAST DB, return summary of matches"""
        fasta = f">query\n{self.sequence}\n"

        with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta") as query_file:
            query_file.write(fasta)
            query_file.flush()

            cmd = [
                "blastp",
                "-query",
                query_file.name,
                "-db",
                str(settings.data_dir.joinpath("blastlib/mite_blastfiles").resolve()),
                "-evalue",
                f"1e-{self.e_val}",
                "-outfmt",
                "5",
            ]

            result = subprocess.run(cmd, capture_output=True, check=True, timeout=5)
            blast_record = Blast.read(io.BytesIO(result.stdout))

            summary = {}
            for hit in blast_record:
                key = hit[0].target.description.split()[0]
                summary[key] = {
                    "sequence_similarity": round(
                        (
                            float(hit[0].annotations.get("positive"))
                            / float(len(self.sequence))
                        )
                        * 100,
                        0,
                    ),
                    "alignment_score": round(hit[0].score, 0),
                    "evalue": round(hit[0].annotations.get("evalue"), 0),
                    "bit_score": round(hit[0].annotations.get("bit score"), 0),
                }

            return summary


# TODO: Create DatabaseQuery class, move all database stuff there


class QueryManager(BaseModel):
    """Organize functions for (database) querying

    Attributes:
        accessions: MITE accession for filtering
        entries: MITE entries
        headers: overview table headers
        operators: translates querybuilder operators to SQLAlchemy col, val pairs
        field_map: safeguard against illegal fields
    """

    accessions: set
    entries: dict
    headers: list
    operators: ClassVar[dict[str]] = {
        "equal": lambda col, val: col == val,
        "not_equal": lambda col, val: col != val,
        "contains": lambda col, val: col.ilike(f"%{val}%"),
        "not_contains": lambda col, val: ~col.ilike(f"%{val}%"),
        "is_null": lambda col, _: col.is_(None),
        "is_not_null": lambda col, _: col.is_not(None),
    }
    field_map: ClassVar[set[str]] = {
        "orcids" "references",
        "evidences",
        "tailoring",
        "enzyme.name",
        "enzyme.enzyme_description",
        "enzyme.uniprot_id",
        "enzyme.genpept_id",
        "enzyme.mibig_id",
        "enzyme.wikidata_id",
        "enzyme.has_auxenzymes",
        "enzyme.organism_id",
        "enzyme.domain_id",
        "enzyme.kingdom_id",
        "enzyme.phylum_id",
        "enzyme.class_id",
        "enzyme.order_id",
        "enzyme.family_id",
        "enzyme.cofactors",
        "reactions.description",
        "reactions.reaction_smarts",
        "reactions.rhea_id",
        "reactions.ec_id",
        "reactions.example_reactions.smiles_substrate",
        "reactions.example_reactions.is_intermediate",
        "reactions.example_reactions.products.smiles_product",
    }

    def stream_csv(self):
        """Convert and stream the (filtered) list of entries"""
        if not self.entries:
            return

        rows = self.get_entries()

        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=rows[0].keys())
        writer.writeheader()

        for row in rows:
            buffer.seek(0)
            buffer.truncate(0)
            writer.writerow(row)
            yield buffer.getvalue()

    def get_entries(self) -> list[dict[str]]:
        """Return (filtered) list of entries"""
        return [v for k, v in self.entries.items() if k in self.accessions]

    def query_db(self, forms: dict, db: Session):
        """Queries the DB

        Args:
            forms: the raw form from QueryBuilder JSON rules
            db: the database session

        Returns:
            A set of MITE accession IDs for filtering
        """
        rules = json.loads(forms["rules"])
        if not rules or len(rules["rules"]) == 0:
            return

        filters = self.parse_rules_to_filters(rules, Entry)
        query = db.query(Entry).filter(filters)
        self.accessions.intersection_update({e.accession for e in query})
        return

    def parse_rules_to_filters(self, rules: dict, base_model: Any) -> ColumnElement:
        """Convert QueryBuilder JSON rules into a SQLAlchemy-compatible expression

        Note: all expressions are chained together with either 'AND' or 'OR'
        """

        expressions = []
        for rule in rules.get("rules", []):
            field_path = rule["field"] if rule["field"] in self.field_map else None
            if not field_path:
                continue

            filter_expr = self.build_filter_from_path(
                model=base_model,
                path_parts=field_path.split("."),
                operator=rule["operator"],
                value=rule.get("value"),
            )
            expressions.append(filter_expr)

        if rules.get("condition", "AND").upper() == "OR":
            return or_(*expressions)
        else:
            return and_(*expressions)

    def build_filter_from_path(
        self, model: Any, path_parts: list, operator: str, value: str
    ):
        """Build a SQLAlchemy filter from a dotted path."""
        mapper = inspect(model)

        if path_parts[0] in mapper.relationships:
            rel = mapper.relationships[path_parts[0]]
            related_model = rel.entity.class_

            nested_filter = self.build_filter_from_path(
                related_model, path_parts[1:], operator, value
            )

            if rel.uselist:
                return rel.class_attribute.any(nested_filter)
            else:
                return rel.class_attribute.has(nested_filter)
        elif path_parts[0] in mapper.columns:
            col = mapper.columns[path_parts[0]]
            return self.operators[operator](col, value)
        else:
            raise ValueError(f"Invalid path segment: {path_parts[0]}")

    def query_sequence(self, forms: dict) -> None:
        if not forms.get("sequence") or not forms.get("e_val"):
            return

        mgr = BlastManager(
            sequence=forms.get("sequence"),
            e_val=forms.get("e_val"),
        )
        results = mgr.run_blastp()
        self.accessions.intersection_update({k for k in results})
        for k, v in results.items():
            self.entries[k].update(v)

        self.headers.extend(
            [
                ("evalue", "E-Value"),
                ("sequence_similarity", "Sequence Sim. (%)"),
                ("alignment_score", "Alignment Score"),
                ("bit_score", "Bit-Score"),
            ]
        )

        return
        # TODO: update headers

        # TODO: continue here
