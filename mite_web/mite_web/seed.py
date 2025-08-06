import json
import os

from flask import current_app

from mite_web.config.extensions import db
from mite_web.models import ChangeLog, Cofactor, Entry, Enzyme, Person, Reference

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
                retirement_reasons=data.get("retirementReasons"),
                comment=data.get("comment"),
            )

            entry.changelogs = get_changelogs(entry, data["changelog"])
            entry.enzyme = get_enzyme(entry, data["enzyme"])

            db.session.add(entry)

    db.session.commit()
    print(f"Seeded database.")


def get_changelogs(entry: Entry, logs: list) -> list[ChangeLog]:
    """Parse changelog and create many-to-many person table"""

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


def get_or_create_referece(doi: str) -> Reference:
    """Add reference if not already existing"""
    reference = Reference.query.filter_by(doi=doi).first()
    if not reference:
        reference = Reference(doi=doi)
        db.session.add(reference)
    return reference


def get_enzyme(entry: Entry, data: dict) -> Enzyme:
    """Parse enzyme info and create many-to-many Cofactor and Reference tables"""

    def _get_or_create_cofactor(cfname: str, cftipo: str) -> Cofactor:
        cofactor = Cofactor.query.filter_by(cofactor_name=cfname).first()
        if not cofactor:
            cofactor = Cofactor(cofactor_name=cfname, cofactor_type=cftipo)
            db.session.add(cofactor)
        return cofactor

    enzyme = Enzyme(
        name=data["name"],
        enzyme_description=data.get("description"),
        uniprot_id=data.get("databaseIds", {}).get("uniprot"),
        genpept_id=data.get("databaseIds", {}).get("genpept"),
        mibig_id=data.get("databaseIds", {}).get("mibig"),
        wikidata_id=data.get("databaseIds", {}).get("wikidata"),
        has_auxenzymes=bool(data.get("auxiliaryEnzymes")),
        entry=entry,
    )
    db.session.add(enzyme)

    enzyme.cofactors = []
    for tipo in ("inorganic", "organic"):
        for name in data.get("cofactors", {}).get(tipo, []):
            enzyme.cofactors.append(_get_or_create_cofactor(name, tipo))

    enzyme.references = []
    enzyme.references.extend(
        [get_or_create_referece(ref) for ref in data["references"]]
    )

    return enzyme
