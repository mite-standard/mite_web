from sqlalchemy import Column, Integer, String, Text

from app.database import Base


class Entries(Base):
    """Stores relational data on baked-in dataset"""

    __tablename__ = "entries"

    accession = Column(String, primary_key=True)
    orcids = Column(Text)
    references = Column(Text)
    evidences = Column(Text)
    tailoring = Column(Text)

    # TODO: Add rest of model
