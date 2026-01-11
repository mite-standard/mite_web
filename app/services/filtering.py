import csv
import io
import json

from pydantic import BaseModel
from sqlalchemy.orm import Session

from app.db.crud import search_db
from app.schemas.blast import BlastSearch
from app.schemas.reaction import ReactionSearch
from app.schemas.structure import StructureSearch


class FilterManager(BaseModel):
    """Organize functions for (database) querying

    Attributes:
        accessions: MITE accession for filtering
        entries: MITE entries
        headers: overview table headers
    """

    accessions: set
    entries: dict
    headers: list

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

    def query_db(self, forms: dict, db: Session) -> None:
        """Query the DB using QueryBuilder query"""
        rules = json.loads(forms["rules"])
        if not rules or len(rules["rules"]) == 0:
            return

        self.accessions.intersection_update(search_db(rules=rules, db=db))

    def query_sequence(self, forms: dict) -> None:
        """Query a sequence against MITE protein sequences"""
        if not forms.get("sequence") or not forms.get("e_val"):
            return

        search = BlastSearch(
            sequence=forms.get("sequence"),
            e_val=forms.get("e_val"),
        )
        results = search.run_blastp()
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

    def query_structure(self, forms: dict) -> None:
        """Query a (sub)structure as smiles/smarts against MITE substrate/product smiles"""
        if not forms.get("substructure_query") or not forms.get("substructure_type"):
            return

        search = StructureSearch(
            structure=forms.get("substructure_query"),
            s_type=forms.get("substructure_type"),
        )
        results = search.search_structure()
        self.accessions.intersection_update({k for k in results})
        for k, v in results.items():
            self.entries[k].update(v)

        self.headers.extend(
            [
                ("reaction", "Reaction"),
                ("example", "Example Reaction"),
            ]
        )
        return

    def query_reaction(self, forms: dict) -> None:
        """Query a reaction smarts against MITE reaction smarts"""
        if (
            not forms.get("reaction_query")
            or not forms.get("similarity")
            or not forms.get("reaction_smarts_radio")
        ):
            return

        search = ReactionSearch(
            reaction=forms.get("reaction_query"),
            sim_score=forms.get("similarity"),
            mode=forms.get("reaction_smarts_radio"),
        )
        results = search.search_reaction()
        self.accessions.intersection_update({k for k in results})
        for k, v in results.items():
            self.entries[k].update(v)

        self.headers.extend(
            [
                ("reaction", "Reaction"),
                ("sim_score", "Similarity Score"),
            ]
        )
        return
