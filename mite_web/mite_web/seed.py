import json
import os

from flask import current_app

from mite_web.config.extensions import db
from mite_web.models import (
    Cofactor,
    Entry,
    Enzyme,
    Evidence,
    ExampleReaction,
    Person,
    Product,
    Reaction,
    Reference,
    Tailoring,
)


def seed_data() -> None:
    """Seeds database with MITE data"""
    if Entry.query.first():
        current_app.logger.warning("Database already seeded.")
        return

    src = current_app.config["DATA_JSON"]
    for entry in src.iterdir():
        with open(entry) as infile:
            data = json.load(infile)

            if data["status"] != "active":
                current_app.logger.warning(
                    f"{data["accession"]} is retired - database parsing skipped."
                )
                continue

            entry = Entry(
                accession=data["accession"],
                orcids=get_orcids(data["changelog"])
            )

            entry.enzyme = get_enzyme(entry, data["enzyme"], current_app.config["SUMMARY"].get(data["accession"], {}))
            entry.reactions = get_reactions(entry, data["reactions"])

            db.session.add(entry)

    db.session.commit()
    current_app.logger.info("Database seeding completed.")


def get_orcids(logs: list) -> str:
    """Parse changelog and concatenate ORCIDs in a string"""
    seen_orcids = set()
    for log in logs:
        for orcid in log.get("contributors", []) + log.get("reviewers", []):
            if orcid not in seen_orcids:
                seen_orcids.add(orcid)
    return "|" + "|".join(seen_orcids) + "|"


def get_or_create_reference(doi: str) -> Reference:
    """Add reference if not already existing"""
    reference = Reference.query.filter_by(doi=doi).first()
    if not reference:
        reference = Reference(doi=doi)
        db.session.add(reference)
    return reference


def get_enzyme(entry: Entry, data: dict, summary: dict) -> Enzyme:
    """Parse enzyme info and create many-to-many Cofactor and Reference tables

    Args:
        entry: an Entry instance
        data: a MITE JSON enzyme dict
        summary: the summary of the MITE entry containing taxonomy info
    """
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
        organism_id=summary.get("organism"),
        domain_id=summary.get("domain"),
        kingdom_id=summary.get("kingdom"),
        phylum_id=summary.get("phylum"),
        class_id=summary.get("class"),
        order_id=summary.get("order"),
        family_id=summary.get("family"),
        entry=entry,
    )
    db.session.add(enzyme)

    enzyme.cofactors = []
    for tipo in ("inorganic", "organic"):
        for name in data.get("cofactors", {}).get(tipo, []):
            enzyme.cofactors.append(_get_or_create_cofactor(name, tipo))

    enzyme.references = []
    enzyme.references.extend(
        [get_or_create_reference(ref) for ref in data["references"]]
    )

    return enzyme


def get_reactions(entry: Entry, data: list) -> list[Reaction]:
    """Parse enzyme info and create many-to-many Cofactor and Reference tables"""

    def _get_or_create_tailoring(item: str) -> Tailoring:
        tailoring = Tailoring.query.filter_by(tailoring=item).first()
        if not tailoring:
            tailoring = Tailoring(tailoring=item)
            db.session.add(tailoring)
        return tailoring

    def _get_or_create_evidence(item: str) -> Evidence:
        evidence = Evidence.query.filter_by(evidence=item).first()
        if not evidence:
            evidence = Evidence(evidence=item)
            db.session.add(evidence)
        return evidence

    def _create_example_reaction(item: dict, r_ref: Reaction) -> ExampleReaction:
        example_reaction = ExampleReaction(
            smiles_substrate=item["substrate"],
            is_intermediate=item["isIntermediate"],
            products=[Product(smiles_product=s) for s in item["products"]],
            reaction=r_ref,
        )
        db.session.add(example_reaction)
        return example_reaction

    reactions = []
    for r in data:
        reaction = Reaction(
            description=r.get("description"),
            reaction_smarts=r["reactionSMARTS"],
            entry=entry,
        )
        db.session.add(reaction)

        reaction.tailorings = [_get_or_create_tailoring(t) for t in r["tailoring"]]
        reaction.evidences = [
            _get_or_create_evidence(e) for e in r["evidence"]["evidenceCode"]
        ]
        reaction.references = []
        reaction.references.extend(
            [get_or_create_reference(ref) for ref in r["evidence"]["references"]]
        )
        reaction.example_reactions = [
            _create_example_reaction(expl, reaction) for expl in r["reactions"]
        ]

    return reactions
