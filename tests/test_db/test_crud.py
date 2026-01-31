from pathlib import Path

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from app.db.crud import enzyme_exists

db_location = Path("tests/dummydata/app.db").resolve()


@pytest.fixture
def db_session():
    engine = create_engine(
        f"sqlite:///{db_location}?mode=ro",
        connect_args={"check_same_thread": False, "uri": True},
    )
    Session = sessionmaker(bind=engine, autoflush=False, autocommit=False)
    session = Session()

    yield session

    session.close()


def test_enzyme_exists_valid(db_session):
    assert enzyme_exists(db=db_session, accession="Q194P6") == "MITE0000020"
    assert enzyme_exists(db=db_session, accession="CAK50792.1") == "MITE0000020"


def test_enzyme_exists_invalid(db_session):
    assert enzyme_exists(db=db_session, accession="asiunasnfasoi") is None
