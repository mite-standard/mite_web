from sqlalchemy.orm import Session

from app.db.models import Entry


def get_entires(db: Session):
    return db.query(Entry).all()
