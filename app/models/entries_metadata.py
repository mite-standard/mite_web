from sqlalchemy import Column, DateTime, Integer, String, Text

from app.database import Base


class EntriesMetadata(Base):
    """Stores metadata on baked-in dataset"""

    __tablename__ = "entries_metadata"

    id = Column(Integer, primary_key=True)
    version = Column(String, nullable=False)
    timestamp = Column(DateTime, nullable=False)
