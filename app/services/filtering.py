import csv
import io
import json

from pydantic import BaseModel
from sqlalchemy.orm import Session

from app.db.crud import query_db
from app.schemas.blast import BlastManager


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

        self.accessions.intersection_update(query_db(rules=rules, db=db))

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
