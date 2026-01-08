from sqlalchemy.orm import Session

from app.models.entries import Entries


def get_entires(db: Session):
    return db.query(Entries).all()
