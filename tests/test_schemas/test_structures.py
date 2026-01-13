import pytest

from app.schemas.structures import GetStructure, SmilesParams


def test_chirality_chiral_input():
    params = SmilesParams(smiles="C[C@H](N)C(=O)O")
    structure = GetStructure()
    structure.strip_chirality(p=params)
    assert structure.smiles == "CC(N)C(=O)O"


def test_chirality_achiral_input():
    params = SmilesParams(smiles="CCCCC")
    structure = GetStructure()
    structure.strip_chirality(p=params)
    assert structure.smiles == "CCCCC"
