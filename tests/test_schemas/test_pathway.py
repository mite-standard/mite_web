import pytest

from app.schemas.pathway import PathwayParams


def test_init_valid():
    params = PathwayParams(substrate="CCCC", mite_id=["MITE0000001"], reaction_nr=[1])
    assert params.reactions == [("MITE0000001", 1)]


def test_init_invalid_smiles():
    with pytest.raises(ValueError):
        PathwayParams(substrate="ertzuiop", mite_id=["MITE0000001"], reaction_nr=[1])


def test_init_invalid_int():
    with pytest.raises(ValueError):
        PathwayParams(substrate="CCCC", mite_id=["MITE0000001"], reaction_nr=[-1])


def test_init_invalid_accession():
    with pytest.raises(ValueError):
        PathwayParams(substrate="CCCC", mite_id=["asdasfdasd"], reaction_nr=[1])


def test_init_invalid_len():
    with pytest.raises(ValueError):
        PathwayParams(substrate="CCCC", mite_id=["asdasfdasd"], reaction_nr=[])
