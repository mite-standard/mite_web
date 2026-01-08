import logging
import sys
from pathlib import Path

from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from app.db import models
from app.db.database import Base

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
if not logger.handlers:
    logger.addHandler(handler)

DATA_DIR = Path(__file__).parent.parent.joinpath("data")
DB_PATH = DATA_DIR / "app.db"

# TODO: Add the full seeding; also add the full models


def main():
    """Create and populate SQLite DB for Docker image."""

    if not DATA_DIR.exists():
        raise RuntimeError(f"Data directory {DATA_DIR} does not exist - Abort.")

    if DB_PATH.exists():
        logger.warning("Existing database found — recreating")

    engine = create_engine(f"sqlite:///{DB_PATH.as_posix()}")

    logger.info("Dropping existing tables")
    Base.metadata.drop_all(engine)

    logger.info("Creating tables")
    Base.metadata.create_all(engine)

    logger.info("Seeding data")
    with Session(engine) as session:
        try:
            session.add(models.Entries(accession="MITE0000002"))
            session.commit()

            # TODO: only for debug - remove this part
            check = session.query(models.Entries).all()
            print({e.accession for e in check})
        except Exception:
            session.rollback()
            raise

    logger.info(f"Databases successfully created at {DB_PATH}")


if __name__ == "__main__":
    main()
