from pathlib import Path

import pytest

from app.services import validation


@pytest.mark.parametrize(
    "v",
    ["0000-0001-6534-6609", "AAAAAAAAAAAAAAAAAAAAAAAA"],
)
def test_submitter_valid(v):
    assert validation.validate_submitter(v) is None


@pytest.mark.parametrize(
    "v",
    ["0000-0001-6534-6609ASA", "asdfasd"],
)
def test_submitter_invalid(v):
    with pytest.raises(ValueError):
        validation.validate_submitter(v)


def test_exists_valid():
    assert validation.check_exists(Path("tests/dummydata/MITE0000001.json")) is None


def test_exists_invalid():
    with pytest.raises(FileNotFoundError):
        validation.check_exists(Path("tests/dummydata/MITE0000000.json"))


@pytest.mark.slow
def test_genpept_valid():
    assert validation.is_genpept(acc="AAK83184.1", timeout=10)


@pytest.mark.slow
def test_genpept_invalid():
    assert validation.is_genpept(acc="asfassad", timeout=10) is False


def test_genpept_timeout():
    assert validation.is_genpept(acc="AAK83184.1", timeout=0.000000001) is None


@pytest.mark.slow
def test_uniprot_valid():
    assert validation.is_uniprot(acc="Q93KW1", timeout=10)


@pytest.mark.slow
def test_uniparc_valid():
    assert validation.is_uniprot(acc="UPI0005831253", timeout=10)


@pytest.mark.slow
def test_uniprot_valid():
    assert validation.is_uniprot(acc="AAK83184.1", timeout=10) is False


def test_uniprot_timeout():
    assert validation.is_uniprot(acc="Q93KW1", timeout=0.000000001) is None


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(container=["1", "2"], msg="message"),
        dict(container=["1"], msg="message"),
    ],
)
def test_nr_items_valid(kwargs):
    assert validation.validate_nr_items(**kwargs) is None


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(container=[], msg="message"),
        dict(container=["", ""], msg="message"),
    ],
)
def test_nr_items_invalid(kwargs):
    with pytest.raises(ValueError):
        validation.validate_nr_items(**kwargs)
