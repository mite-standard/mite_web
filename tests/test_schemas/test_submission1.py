from pathlib import Path

import pytest
from pydantic_core._pydantic_core import ValidationError

from app.schemas.submission import (
    ExistDraftData,
    ExistDraftForm,
    ExistDraftService,
    NewDraftData,
    NewDraftForm,
    NewDraftService,
)


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
