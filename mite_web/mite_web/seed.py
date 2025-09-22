import json
import os
import re

from flask import current_app

from mite_web.config.extensions import db
from mite_web.models import (
    Entry,
    Enzyme,
    ExampleReaction,
    Product,
    Reaction,
)


def seed_data() -> None:
    """Seeds database with MITE data"""
    if Entry.query.first():
        current_app.logger.warning("Database already seeded.")
        return

    src = current_app.config["DATA_JSON"]
    for entry in src.iterdir():
        if not re.fullmatch(r"MITE[0-9]{7}.json", entry.name):
            continue

        with open(entry) as infile:
            data = json.load(infile)

            if data["status"] != "active":
                current_app.logger.warning(
                    f"{data["accession"]} is retired - database parsing skipped."
                )
                continue

            entry = Entry(
                accession=data["accession"],
                orcids=get_orcids(data["changelog"]),
                references=get_references(data),
                evidences=get_evidence(data["reactions"]),
                tailoring=get_tailoring(data["reactions"]),
            )

            entry.enzyme = get_enzyme(
                entry,
                data["enzyme"],
                current_app.config["SUMMARY_ACTIVE"].get(data["accession"], {}),
            )
            entry.reactions = get_reactions(entry, data["reactions"])

            db.session.add(entry)

    db.session.commit()
    current_app.logger.info("Database seeding completed.")


def get_orcids(logs: list) -> str:
    """Parse changelog and concatenate ORCIDs in a string"""
    seen_orcids = set()
    for log in logs:
        for orcid in log.get("contributors", []) + log.get("reviewers", []):
            seen_orcids.add(orcid)
    return "|".join(seen_orcids)


def get_references(data: dict) -> str:
    """Parse data and concatenate references in a string"""
    seen_refs = set()
    for ref in data["enzyme"]["references"]:
        seen_refs.add(ref.strip("doi:"))

    for r in data["reactions"]:
        for ref in r["evidence"]["references"]:
            seen_refs.add(ref.strip("doi:"))

    return "|".join(seen_refs)


def get_evidence(reactions: list) -> str:
    """Parse data and concatenate evidence"""
    seen_evidence = set()
    for r in reactions:
        for evidence in r["evidence"]["evidenceCode"]:
            seen_evidence.add(evidence)
    return "|".join(seen_evidence)


def get_tailoring(reactions: list) -> str:
    """Parse data and concatenate tailoring terms"""
    seen_tailoring = set()
    for r in reactions:
        for tailoring in r["tailoring"]:
            seen_tailoring.add(tailoring)
    return "|".join(seen_tailoring)


def get_enzyme(entry: Entry, data: dict, summary: dict) -> Enzyme:
    """Parse enzyme info and create many-to-many Cofactor and Reference tables

    Args:
        entry: an Entry instance
        data: a MITE JSON enzyme dict
        summary: the summary of the MITE entry containing taxonomy info
    """

    def _add_cofactors() -> str:
        """Returns concatenated cofactors"""
        cofactor_set = set()
        for cat in ("inorganic", "organic"):
            for name in data.get("cofactors", {}).get(cat, []):
                cofactor_set.add(name)
        return "|".join(cofactor_set)

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
        cofactors=_add_cofactors(),
        entry=entry,
    )
    db.session.add(enzyme)

    return enzyme


def get_reactions(entry: Entry, data: list) -> list[Reaction]:
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
            rhea_id=r.get("databaseIds", {}).get("rhea"),
            ec_id=r.get("databaseIds", {}).get("ec"),
            entry=entry,
        )
        db.session.add(reaction)

        reaction.example_reactions = [
            _create_example_reaction(expl, reaction) for expl in r["reactions"]
        ]
        reactions.append(reaction)

    return reactions
