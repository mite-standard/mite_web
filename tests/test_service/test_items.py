from pathlib import Path

import pytest

from app.services.items import MiteModel


@pytest.fixture
def mite_model(monkeypatch) -> MiteModel:
    monkeypatch.setattr(MiteModel, "data_dir", Path("tests/dummydata"))
    return MiteModel(mite_id="MITE0000001")


def test_mite_number(mite_model):
    assert mite_model.mite_number == 1


def test_previous_id(mite_model):
    assert mite_model.previous_id is None


def test_next_id(mite_model):
    assert mite_model.next_id == "MITE0000002"


def test_previous_exists(mite_model):
    assert mite_model.previous_exists() == False


def test_next_exists(mite_model):
    assert mite_model.next_exists() == False


def test_model_invalid_id():
    with pytest.raises(ValueError):
        MiteModel(mite_id="dfhjk")


def test_model_invalid_file(monkeypatch):
    monkeypatch.setattr(MiteModel, "data_dir", Path("tests/dummydata"))
    with pytest.raises(FileNotFoundError):
        MiteModel(mite_id="MITE0000002")
