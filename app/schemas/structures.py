import re

from pydantic import BaseModel, Field, field_validator
from rdkit import Chem


class PeptideParams(BaseModel):
    sequence: str = Field(..., min_length=1)

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


class SmilesParams(BaseModel):
    smiles: str = Field(..., min_length=1)

    @field_validator("smiles")
    @classmethod
    def validate_smiles(cls, v: str):
        if not Chem.MolFromSmiles(v):
            raise ValueError("Input SMILES string invalid according to RDKit.")
        return v


class GetStructure(BaseModel):
    smiles: str = None

    def get_peptide(self, p: PeptideParams):
        """Create smiles from peptide sequence"""
        self.smiles = Chem.MolToSmiles(Chem.MolFromSequence(p.sequence))

    def get_canonical(self, p: SmilesParams):
        """Canonicalize a smiles sequence"""
        self.smiles = Chem.MolToSmiles(Chem.MolFromSmiles(p.smiles))

    def strip_chirality(self, p: SmilesParams):
        """Strip stereochemical information from SMILES string"""
        mol = Chem.MolFromSmiles(p.smiles)
        Chem.RemoveStereochemistry(mol)
        self.smiles = Chem.MolToSmiles(mol)
