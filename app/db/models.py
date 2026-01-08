from sqlalchemy import Column, Integer, String, Text

from app.db.database import Base


class Entries(Base):
    """Stores relational data on baked-in dataset"""

    __tablename__ = "entries"

    accession = Column(String, primary_key=True)

    # TODO: Add rest of model
