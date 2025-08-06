import json
import os

from flask import current_app

from mite_web.config.extensions import db
from mite_web.models import ChangeLog, Entry, Person

# TODO(MMZ 5.8): replace print with logging


def seed_data() -> None:
    """Seeds database with MITE data"""
    if Entry.query.first():
        print("Database already seeded.")
        return

    src = current_app.config["DATA_JSON"]
    for entry in src.iterdir():
        with open(entry) as infile:
            data = json.load(infile)

            entry = Entry(
                accession=data["accession"],
                status=data["status"],
                retirement_reasons=data.get("retirementReasons", ""),
                comment=data.get("comment", ""),
            )

            entry.changelogs = get_changelogs(entry, data["changelog"])

            db.session.add(entry)

    db.session.commit()
    print(f"Seeded databbase.")


def get_changelogs(entry: Entry, logs: list) -> list[ChangeLog]:
    def _get_or_create_person(orcid: str) -> Person:
        person = Person.query.filter_by(orcid=orcid).first()
        if not person:
            person = Person(orcid=orcid)
            db.session.add(person)
        return person

    changelogs = []
    for i in logs:
        cl = ChangeLog(
            version=i["version"], date=i["date"], comment=i["comment"], entry=entry
        )
        db.session.add(cl)

        cl.contributors = [_get_or_create_person(orcid) for orcid in i["contributors"]]
        cl.reviewers = [_get_or_create_person(orcid) for orcid in i["reviewers"]]
        changelogs.append(cl)

    return changelogs
