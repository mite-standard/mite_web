import io
import re
import subprocess
import tempfile

from Bio import Blast
from pydantic import BaseModel, Field, field_validator

from app.core.config import settings


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
