import re
from pathlib import Path

import pytest
from pydantic_core._pydantic_core import ValidationError

from app.schemas.submission import (
    ExistDraftData,
    ExistDraftForm,
    ExistDraftService,
    MiteService,
    NewDraftData,
    NewDraftForm,
    NewDraftService,
)


@pytest.fixture
def form_input():
    return {
        "reaction[0]smarts": [""],
        "reaction[0]description": [""],
        "reaction[0]rhea": [""],
        "reaction[0]ec": [""],
        "reaction[0]tailoring[]": [""],
        "reaction[0]evidencecode[]": [""],
        "reaction[0]ref[]": [""],
        "reaction[0]knownreaction[0]substrate": [""],
        "reaction[0]knownreaction[0]description": [""],
        "reaction[0]knownreaction[0]intermediate": [""],
        "reaction[0]knownreaction[0]products[]": [""],
        "reaction[1]smarts": [""],
        "reaction[1]description": [""],
        "reaction[1]rhea": [""],
        "reaction[1]ec": [""],
        "reaction[1]tailoring[]": [""],
        "reaction[1]evidencecode[]": [""],
        "reaction[1]ref[]": [""],
        "reaction[1]knownreaction[0]substrate": [""],
        "reaction[1]knownreaction[0]description": [""],
        "reaction[1]knownreaction[0]intermediate": [""],
        "reaction[1]knownreaction[0]products[]": [""],
    }


@pytest.mark.slow
@pytest.mark.parametrize(
    "kwargs",
    [
        dict(
            name="AbcD",
            submitter="0000-0001-6534-6609",
            accession="Q93KW1",
            db="uniprot",
            token="str",
        ),
        dict(
            name="AbcD",
            submitter="AAAAAAAAAAAAAAAAAAAAAAAA",
            token="str",
            accession="AAK83184.1",
            db="genpept",
        ),
    ],
)
def test_draft_init_valid(kwargs):
    form = NewDraftForm(**kwargs)
    assert isinstance(form, NewDraftForm)


@pytest.mark.slow
@pytest.mark.parametrize(
    "kwargs",
    [
        dict(),
        dict(
            name="AbcD",
            submitter="https://orcid.org/0000-0001-6534-6609",
            token="str",
        ),
        dict(name="", submitter="", token="str"),
        dict(
            name="AbcD",
            submitter="AAAAAAAAAAAAAAAAAAAAAAAA",
            token="str",
            accession="AAK83184.1",
            db="uniprot",
        ),
    ],
)
def test_draft_init_invalid(kwargs):
    with pytest.raises(ValidationError):
        NewDraftForm(**kwargs)


@pytest.mark.slow
@pytest.mark.parametrize(
    "kwargs",
    [
        dict(
            name="AbcD",
            submitter="0000-0001-6534-6609",
            accession="Q93KW1",
            db="uniprot",
            token="str",
        )
    ],
)
def test_draft_service_valid(kwargs):
    service = NewDraftService()
    results = service.parse(NewDraftForm(**kwargs))
    assert isinstance(results, NewDraftData)
    assert results.data["enzyme"]["databaseIds"]["uniprot"] == "Q93KW1"


@pytest.mark.slow
@pytest.mark.parametrize(
    "kwargs",
    [
        dict(
            name="AbcD",
            submitter="0000-0001-6534-6609",
            accession="safasdfasdfasd",
            db="uniprot",
            token="str",
        )
    ],
)
def test_draft_service_invalid(kwargs):
    with pytest.raises(ValueError):
        service = NewDraftService()
        service.parse(NewDraftForm(**kwargs))


@pytest.fixture
def exist_model(monkeypatch):
    monkeypatch.setattr(ExistDraftForm, "data_dir", Path("tests/dummydata"))


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(
            mite_id="MITE0000001",
            submitter="0000-0001-6534-6609",
            token="str",
        ),
        dict(
            mite_id="MITE0000001",
            submitter="AAAAAAAAAAAAAAAAAAAAAAAA",
            token="str",
        ),
    ],
)
def test_exists_init_valid(exist_model, kwargs):
    form = ExistDraftForm(**kwargs)
    assert isinstance(form, ExistDraftForm)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(
            mite_id="adasdoads",
            submitter="0000-0001-6534-6609",
            token="str",
        )
    ],
)
def test_exists_init_invalid(exist_model, kwargs):
    with pytest.raises(ValidationError):
        ExistDraftForm(**kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(
            mite_id="MITE0000000",
            submitter="AAAAAAAAAAAAAAAAAAAAAAAA",
            token="str",
        ),
    ],
)
def test_exists_init_invalid(exist_model, kwargs):
    with pytest.raises(FileNotFoundError):
        ExistDraftForm(**kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(
            mite_id="MITE0000001",
            submitter="0000-0001-6534-6609",
            token="str",
        ),
    ],
)
def test_exists_service_valid(exist_model, kwargs):
    service = ExistDraftService()
    results = service.parse(ExistDraftForm(**kwargs))
    assert isinstance(results, ExistDraftData)
    assert len(results.data["changelog"]) == 3


def test_parse_reactions(form_input):
    model = MiteService()
    parsed = model.reactions(form_input)
    assert len(parsed) == 2


def test_get_form_instance(form_input):
    model = MiteService()
    parsed = model.get_form_instance(
        form_input, pattern=re.compile(r"reaction\[(\d+)\]")
    )
    assert parsed == [0, 1]


def test_get_known_reactions(form_input):
    model = MiteService()
    parsed = model.get_instance_known_reactions(form_input, idx=0)
    assert parsed == [0]


def test_get_krxs(form_input):
    model = MiteService()
    parsed = model.get_krxs(form_input, rx=0, krx=[0])
    assert len(parsed) == 1
