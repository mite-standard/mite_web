import json

import pytest
from fastapi import HTTPException

from app.services.submission import (
    add_reviewer,
    check_own_entry,
    check_schema,
    pop_dummy_reviewer,
    set_active,
)


@pytest.fixture
def mite_entry():
    with open("tests/dummydata/MITE0000001.json") as file_in:
        return json.load(file_in)


def test_check_schema_valid(mite_entry):
    check_schema(mite_entry)


def test_check_schema_invalid():
    with pytest.raises(ValueError):
        check_schema({})


def test_check_own_entry_valid(mite_entry):
    check_own_entry(user="ssdssa", data=mite_entry)


def test_check_own_entry_invalid(mite_entry):
    with pytest.raises(HTTPException):
        check_own_entry(user="0000-0001-6534-6609", data=mite_entry)


def test_pop_dummy_reviewer_valid():
    data = {"changelog": [{"reviewers": ["BBBBBBBBBBBBBBBBBBBBBBBB"]}]}
    data = pop_dummy_reviewer(data)
    assert "BBBBBBBBBBBBBBBBBBBBBBBB" not in data["changelog"][-1]["reviewers"]


def test_add_reviewer():
    data = {"changelog": [{"reviewers": ["abc"]}]}
    data = add_reviewer(data=data, reviewer="xyz")
    assert "xyz" in data["changelog"][-1]["reviewers"]


def test_set_active_valid():
    data = {"status": "pending"}
    data = set_active(data)
    assert data["status"] == "active"


def test_set_active_invalid():
    data = {"status": "retired"}
    data = set_active(data)
    assert data["status"] == "retired"
