import pytest
from pydantic_core._pydantic_core import ValidationError

from app.schemas.submission import NewDraftForm


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(
            name="AbcD",
            submitter="0000-0001-6534-6609",
        ),
        dict(
            name="AbcD",
            submitter="AAAAAAAAAAAAAAAAAAAAAAAA",
        ),
    ],
)
def test_draft_init_valid(kwargs):
    form = NewDraftForm(**kwargs)
    assert isinstance(form, NewDraftForm)


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(),
        dict(
            name="AbcD",
            submitter="https://orcid.org/0000-0001-6534-6609",
        ),
        dict(name="", submitter=""),
    ],
)
def test_draft_init_invalid(kwargs):
    with pytest.raises(ValidationError):
        NewDraftForm(**kwargs)
