import json
import logging
import re
import sys
from pathlib import Path

from pydantic import BaseModel, Field, model_validator
from sqlalchemy import create_engine
from sqlalchemy.orm import Session

from app.db.database import Base
from app.db.models import Entry, Enzyme, ExampleReaction, Product, Reaction

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
if not logger.handlers:
    logger.addHandler(handler)

# TODO: Add the full seeding; also add the full models


class DBSeeder(BaseModel):
    """Prepares database seeding"""

    data_dir: Path = Field(
        default_factory=lambda: Path(__file__).parent.parent.joinpath("data")
    )
    db_path: Path = Field(
        default_factory=lambda: Path(__file__).parent.parent.joinpath("data/app.db")
    )
    data_src: Path = Field(
        default_factory=lambda: Path(__file__).parent.parent.joinpath("data/data")
    )

    @model_validator(mode="after")
    def validate_data(self):
        if not self.data_dir.exists():
            raise RuntimeError(
                f"Data directory {self.data_dir} does not exist - Abort."
            )
        return self

    def prepare_data(self) -> list:
        """Prepares data for database seeding"""

        seed = []
        for item in self.data_src.iterdir():
            if not re.fullmatch(r"MITE[0-9]{7}.json", item.name):
                continue

            with open(item) as infile:
                data = json.load(infile)

                if data["status"] != "active":
                    logger.warning(f"{data["accession"]} is retired - skipped.")
                    continue

                entry = Entry(
                    accession=data["accession"],
                    orcids=self.get_orcids(data["changelog"]),
                    references=self.get_references(data),
                    evidences=self.get_evidence(data["reactions"]),
                    tailoring=self.get_tailoring(data["reactions"]),
                )

                entry.enzyme = self.get_enzyme(
                    entry, data["enzyme"], active.get(data["accession"], {})
                )
                entry.reactions = self.get_reactions(entry, data["reactions"])
                seed.append(entry)

        return seed

    @staticmethod
    def get_orcids(logs: list) -> str:
        """Parse changelog and concatenate ORCIDs in a string"""
        seen_orcids = set()
        for log in logs:
            for orcid in log.get("contributors", []) + log.get("reviewers", []):
                seen_orcids.add(orcid)
        return "|".join(seen_orcids)

    @staticmethod
    def get_references(data: dict) -> str:
        """Parse data and concatenate references in a string"""
        seen_refs = set()
        for ref in data["enzyme"]["references"]:
            seen_refs.add(ref.strip("doi:"))

        for r in data["reactions"]:
            for ref in r["evidence"]["references"]:
                seen_refs.add(ref.strip("doi:"))

        return "|".join(seen_refs)

    @staticmethod
    def get_evidence(reactions: list) -> str:
        """Parse data and concatenate evidence"""
        seen_evidence = set()
        for r in reactions:
            for evidence in r["evidence"]["evidenceCode"]:
                seen_evidence.add(evidence)
        return "|".join(seen_evidence)

    @staticmethod
    def get_tailoring(reactions: list) -> str:
        """Parse data and concatenate tailoring terms"""
        seen_tailoring = set()
        for r in reactions:
            for tailoring in r["tailoring"]:
                seen_tailoring.add(tailoring)
        return "|".join(seen_tailoring)

    @staticmethod
    def get_enzyme(entry: Entry, data: dict, summary: dict, session: Session) -> Enzyme:
        """Parse enzyme info"""

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
        session.add(enzyme)

        return enzyme

    @staticmethod
    def get_reactions(entry: Entry, data: list, session: Session) -> list[Reaction]:
        """Parse reaction info"""

        def _create_example_reaction(item: dict, r_ref: Reaction) -> ExampleReaction:
            example_reaction = ExampleReaction(
                smiles_substrate=item["substrate"],
                is_intermediate=item["isIntermediate"],
                products=[Product(smiles_product=s) for s in item["products"]],
                reaction=r_ref,
            )
            session.add(example_reaction)
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
            session.add(reaction)

            reaction.example_reactions = [
                _create_example_reaction(expl, reaction) for expl in r["reactions"]
            ]
            reactions.append(reaction)

        return reactions


def main():
    """Create and populate SQLite DB for Docker image."""
    seeder = DBSeeder()

    engine = create_engine(f"sqlite:///{seeder.db_path.as_posix()}")

    logger.info("Dropping existing tables")
    Base.metadata.drop_all(engine)

    logger.info("Creating tables")
    Base.metadata.create_all(engine)

    logger.info("Seeding data")

    with open(seeder.data_dir.joinpath("active.json")) as infile:
        active = json.load(infile)
    with Session(engine) as session:
        try:
            for key in active:
                with open(seeder.data_src.joinpath(f"{key}.json")) as infile:
                    data = json.load(infile)

                entry = Entry(
                    accession=data["accession"],
                    orcids=seeder.get_orcids(data["changelog"]),
                    references=seeder.get_references(data),
                    evidences=seeder.get_evidence(data["reactions"]),
                    tailoring=seeder.get_tailoring(data["reactions"]),
                )

                entry.enzyme = seeder.get_enzyme(
                    entry, data["enzyme"], active.get(data["accession"], {}), session
                )
                entry.reactions = seeder.get_reactions(
                    entry, data["reactions"], session
                )

                session.add(entry)

            session.commit()
        except Exception:
            session.rollback()
            raise

    logger.info(f"Databases successfully created at {seeder.db_path}")


if __name__ == "__main__":
    main()
