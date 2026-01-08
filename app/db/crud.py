from sqlalchemy.orm import Session

from app.db.models import Entries


def get_entires(db: Session):
    return db.query(Entries).all()
