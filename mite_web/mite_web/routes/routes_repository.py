"""Routes for repository pages.

Copyright (c) 2024-present Mitja Maximilian Zdouc, PhD

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import copy
import json
import os
import pickle
import re
import subprocess
import uuid
from pathlib import Path
from typing import Any, ClassVar

import pandas as pd
from Bio import Blast, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from flask import current_app, flash, render_template, request
from pydantic import BaseModel
from rdkit.Chem import MolFromSmarts, MolFromSmiles, PandasTools, rdChemReactions
from rdkit.DataStructs import FingerprintSimilarity
from sqlalchemy import and_, inspect, or_

from mite_web.config.extensions import db
from mite_web.models import Entry
from mite_web.routes import bp


class StructureManager(BaseModel):
    """Organize querying functions

    Attributes:
        summary: a dict of the entries
        dump_smiles: path to csv-file of all SMILES in MITE entries
        substrates_pickle: path to pre-calculated fingerprints of substrate SMILES
        products_pickle: path to pre-calculated fingerprints of products SMILES
        reaction_pickle: path to pre-calculated fingerprints of reaction SMARTS
    """

    summary: dict
    dump_smiles: Path = Path(__file__).parent.parent.joinpath(
        "data/download/dump_smiles.csv"
    )
    substrates_pickle: Path = Path(__file__).parent.parent.joinpath(
        "data/substrate_list.pickle"
    )
    products_pickle: Path = Path(__file__).parent.parent.joinpath(
        "data/product_list.pickle"
    )
    reaction_pickle: Path = Path(__file__).parent.parent.joinpath(
        "data/reaction_fps.pickle"
    )

    def return_summary(self) -> dict:
        return self.summary

    def query_substructure(self, action: str, query: str) -> set:
        """Query dataset for a specific substructure and filter summary for matching entries

        Arguments:
            action: the type of matching to perform (smiles or smarts matching)
            query: a SMILES or SMARTS string

        Returns:
            A set of relevant MITE accessions
        """
        with open(
            self.substrates_pickle,
            "rb",
        ) as infile:
            substrate_list = pickle.load(infile)

        with open(self.products_pickle, "rb") as infile:
            product_list = pickle.load(infile)

        df = pd.read_csv(self.dump_smiles)
        df["ROMol_substrates"] = substrate_list
        df["ROMol_products"] = product_list

        if action == "smiles":
            mol_query = MolFromSmiles(query)
        else:
            mol_query = MolFromSmarts(query)

        df_substrate_match = df[df["ROMol_substrates"] >= mol_query]
        df_product_match = df[df["ROMol_products"] >= mol_query]

        matches = set()
        matches.update(set(df_substrate_match["mite_id"]))
        matches.update(set(df_product_match["mite_id"]))

        matches = sorted(matches)

        mite_accessions = set()
        for match in matches:
            key = match.split(".")[0]
            mite_accessions.add(key)
            reaction = match.split(".")[1].replace("reaction", "rx")
            example = match.split(".")[2].replace("example", "ex")

            if self.summary[key].get("reaction"):
                self.summary[key]["reaction"].append(reaction)
                self.summary[key]["example"].append(f"{reaction}-{example}")
            else:
                self.summary[key]["reaction"] = [reaction]
                self.summary[key]["example"] = [f"{reaction}-{example}"]

        return mite_accessions

    def query_reaction(self, query: str, sim_score: float) -> set:
        """Query dataset for a reaction and filter summary for matching entries

        Arguments:
            query: a reaction SMARTS string
            sim_score: the minimum similarity score to consider an entry

        Returns:
            A set of relevant MITE accessions
        """

        # TODO: implement the non-greedy fingerprint
        try:
            sim_score = float(sim_score)
        except Exception as e:
            raise ValueError("Similarity score not a number") from e

        with open(self.reaction_pickle, "rb") as infile:
            reaction_dict = pickle.load(infile)

        query_fp = rdChemReactions.CreateStructuralFingerprintForReaction(
            rdChemReactions.ReactionFromSmarts(query)
        )

        similarities = [
            FingerprintSimilarity(query_fp, lib_fp)
            for lib_fp in reaction_dict["reaction_fps"]
        ]
        reaction_dict["similarities"] = similarities

        mite_accessions = set()
        df = pd.DataFrame(reaction_dict)
        for _, row in df.iterrows():
            key = row["mite_id"].split(".")[0]
            if row["similarities"] >= sim_score:
                mite_accessions.add(key)
                if self.summary[key].get("reaction"):
                    self.summary[key]["reaction"].append(
                        row["mite_id"].split(".")[1].replace("reaction", "rx")
                    )
                    self.summary[key]["sim_score"].append(round(row["similarities"], 2))
                else:
                    self.summary[key]["reaction"] = [
                        row["mite_id"].split(".")[1].replace("reaction", "rx")
                    ]
                    self.summary[key]["sim_score"] = [round(row["similarities"], 2)]

        return mite_accessions


class DatabaseManager:
    """Organize functions for database querying

    Attributes:
        operators: translates querybuilder operators to SQLAlchemy col, val pairs
    """

    operators: ClassVar = {
        "equal": lambda col, val: col == val,
        "not_equal": lambda col, val: col != val,
        "contains": lambda col, val: col.ilike(f"%{val}%"),
        "not_contains": lambda col, val: ~col.ilike(f"%{val}%"),
        "is_null": lambda col, _: col.is_(None),
        "is_not_null": lambda col, _: col.is_not(None),
    }

    field_map: ClassVar = {
        "orcids": "orcids",
        "references": "references",
        "evidences": "evidences",
        "tailoring": "tailoring",
        "enzyme.name": "enzyme.name",
        "enzyme.enzyme_description": "enzyme.enzyme_description",
        "enzyme.uniprot_id": "enzyme.uniprot_id",
        "enzyme.genpept_id": "enzyme.genpept_id",
        "enzyme.mibig_id": "enzyme.mibig_id",
        "enzyme.wikidata_id": "enzyme.wikidata_id",
        "enzyme.has_auxenzymes": "enzyme.has_auxenzymes",
        "enzyme.organism_id": "enzyme.organism_id",
        "enzyme.domain_id": "enzyme.domain_id",
        "enzyme.kingdom_id": "enzyme.kingdom_id",
        "enzyme.phylum_id": "enzyme.phylum_id",
        "enzyme.class_id": "enzyme.class_id",
        "enzyme.order_id": "enzyme.order_id",
        "enzyme.family_id": "enzyme.family_id",
        "enzyme.cofactors": "enzyme.cofactors",
        "reactions.description": "reactions.description",
        "reactions.reaction_smarts": "reactions.reaction_smarts",
        "reactions.rhea_id": "reactions.rhea_id",
        "reactions.ec_id": "reactions.ec_id",
        "reactions.example_reactions.smiles_substrate": "reactions.example_reactions.smiles_substrate",
        "reactions.example_reactions.is_intermediate": "reactions.example_reactions.is_intermediate",
        "reactions.example_reactions.products.smiles_product": "reactions.example_reactions.products.smiles_product",
    }

    def query_db(self, rules: dict) -> set:
        """Queries the DB

        Args:
            rules: the raw query from QueryBuilder JSON rules

        Returns:
            A set of MITE accession IDs for filtering
        """
        filters = self.parse_rules_to_filters(rules, Entry)
        query = db.session.query(Entry).filter(filters)
        return {e.accession for e in query}

    def parse_rules_to_filters(self, rules: dict, base_model: Any) -> list:
        """Convert QueryBuilder JSON rules into a SQLAlchemy-compatible expression"""
        expressions = []

        for rule in rules.get("rules", []):
            field_path = self.field_map.get(rule["field"])
            if not field_path:
                continue

            operator = rule["operator"]
            value = rule.get("value")

            path_parts = field_path.split(".")
            filter_expr = self.build_filter_from_path(
                base_model, path_parts, operator, value
            )
            expressions.append(filter_expr)

        if rules.get("condition", "AND").upper() == "OR":
            return or_(*expressions)
        else:
            return and_(*expressions)

    def build_filter_from_path(
        self, model: Any, path_parts: list, operator: str, value: str
    ):
        """Build a SQLAlchemy filter from a dotted path."""
        mapper = inspect(model)

        if path_parts[0] in mapper.relationships:
            rel = mapper.relationships[path_parts[0]]
            related_model = rel.entity.class_

            nested_filter = self.build_filter_from_path(
                related_model, path_parts[1:], operator, value
            )

            if rel.uselist:
                return rel.class_attribute.any(nested_filter)
            else:
                return rel.class_attribute.has(nested_filter)
        elif path_parts[0] in mapper.columns:
            col = mapper.columns[path_parts[0]]
            return self.operators[operator](col, value)
        else:
            raise ValueError(f"Invalid path segment: {path_parts[0]}")


class BlastManager(BaseModel):
    """Organizes querying for a protein sequence

    Attributes:
        summary: entry summary for overview page
        blastlib: Path to blast library

    """

    summary: dict
    blastlib: Path = Path(__file__).parent.parent.joinpath("data/blastlib/")

    def return_summary(self) -> dict:
        return self.summary

    def query_sequence(self, query: str, e_val: int) -> set:
        """Run BLAST against MITE BLAST DB for a specific protein and filter summary for matching entries

        Arguments:
            query: a protein sequence
            e_val: the e-value to filter matches with

        Returns:
            A set of relevant MITE accessions

        Raises:
            RuntimeError: input sanitization detected illegal input
            ValueError: E-value invalid
        """
        query = query.replace("\n", "")
        query = re.sub(r"\s+", "", query, flags=re.UNICODE)

        try:
            e_val = int(e_val)
        except Exception as e:
            raise ValueError("E-value not an integer") from e

        if not query.isalpha():
            raise RuntimeError(
                "The query AA sequence contains illegal characters - please check."
            )
        else:
            query = query.upper()

        job_uuid = uuid.uuid1()

        record = SeqRecord(seq=Seq(query), id="user_query", description="user_query")

        with open(self.blastlib.joinpath(f"{job_uuid}.fasta"), "w") as outfile:
            SeqIO.write(record, outfile, "fasta")

        command = [
            "blastp",
            "-query",
            f"{self.blastlib.joinpath(f"{job_uuid}.fasta")}",
            "-db",
            f"{self.blastlib.joinpath("mite_blastfiles")}",
            "-out",
            f"{self.blastlib.joinpath(f"{job_uuid}.xml")}",
            "-evalue",
            f"1e-{e_val}",
            "-outfmt",
            "5",
        ]
        subprocess.run(command, check=True, timeout=5)

        with open(self.blastlib.joinpath(f"{job_uuid}.xml"), "rb") as infile:
            blast_record = Blast.read(infile)

        filtered_accessions = set()
        for hit in blast_record:
            key = hit[0].target.description.split()[0]
            filtered_accessions.add(key)
            self.summary[key]["sequence_similarity"] = round(
                (float(hit[0].annotations.get("positive")) / float(len(query))) * 100, 0
            )
            self.summary[key]["alignment_score"] = round(hit[0].score, 0)
            self.summary[key]["evalue"] = round(hit[0].annotations.get("evalue"), 0)
            self.summary[key]["bit_score"] = round(
                hit[0].annotations.get("bit score"), 0
            )

        os.remove(self.blastlib.joinpath(f"{job_uuid}.xml"))
        os.remove(self.blastlib.joinpath(f"{job_uuid}.fasta"))

        return filtered_accessions


@bp.route("/overview/", methods=["GET", "POST"])
def overview() -> str:
    """Render the repository overview page of mite_web

    If search null or fails, return the full list

    Returns:
        The overview.html page as string.
    """
    summary = copy.deepcopy(current_app.config["SUMMARY"])
    accessions = copy.deepcopy(current_app.config["ACCESSIONS"])
    headers = [
        ("accession", "MITE Accession"),
        ("status", "Status"),
        ("name", "Enzyme Name"),
        ("tailoring", "Tailoring Reaction"),
        ("description", "Enzyme Description"),
        ("organism", "Organism"),
        ("family", "Family"),
        ("reaction_description", "Reaction Description"),
    ]
    filtered = False

    try:
        if request.method == "POST":
            forms = request.form.to_dict()

            rules = json.loads(forms["rules"])
            if len(rules["rules"]) > 0:
                db_manager = DatabaseManager()
                accessions.intersection_update(db_manager.query_db(rules))

            if forms.get("sequence_query", "") != "":
                blast_manager = BlastManager(summary=summary)
                accessions.intersection_update(
                    blast_manager.query_sequence(
                        query=forms.get("sequence_query"),
                        e_val=forms.get("sequence_similarity"),
                    )
                )
                summary = blast_manager.return_summary()
                headers.extend(
                    [
                        ("evalue", "E-Value"),
                        ("sequence_similarity", "Sequence Sim. (%)"),
                        ("alignment_score", "Alignment Score"),
                        ("bit_score", "Bit-Score"),
                    ]
                )

            if forms.get("substructure_query", "") != "":
                structure_manager = StructureManager(summary=summary)
                accessions.intersection_update(
                    structure_manager.query_substructure(
                        action=forms.get("query_substructure"),
                        query=forms.get("substructure_query"),
                    )
                )
                summary = structure_manager.return_summary()
                headers.extend(
                    [
                        ("reaction", "Reaction"),
                        ("example", "Example Reaction"),
                    ]
                )

            if forms.get("reaction_query", "") != "":
                structure_manager = StructureManager(summary=summary)
                accessions.intersection_update(
                    structure_manager.query_reaction(
                        query=forms.get("reaction_query"),
                        sim_score=forms.get("similarity"),
                    )
                )
                summary = structure_manager.return_summary()
                headers.extend(
                    [
                        ("reaction", "Reaction"),
                        ("sim_score", "Similarity Score"),
                    ]
                )

            filtered = True
            job_id = uuid.uuid1()

            summary = [val for key, val in summary.items() if key in accessions]
            df = pd.DataFrame(summary)
            df.to_csv(
                current_app.config["QUERIES"].joinpath(f"{job_id}.csv"), index=False
            )

            return render_template(
                "overview.html",
                entries=summary,
                headers=headers,
                form_vals=current_app.config["FORM_VALS"],
                filtered=filtered,
                job_id=job_id,
            )
    except Exception as e:
        flash(f"An error occurred during search: '{e!s}'")
        summary = copy.deepcopy(current_app.config["SUMMARY"])
        accessions = copy.deepcopy(current_app.config["ACCESSIONS"])

    summary = [val for key, val in summary.items() if key in accessions]
    return render_template(
        "overview.html",
        entries=summary,
        headers=headers,
        form_vals=current_app.config["FORM_VALS"],
        filtered=filtered,
    )


@bp.route("/repository/<mite_acc>")
def repository(mite_acc: str) -> str:
    """Render the individual pages

    Arguments:
        mite_acc: the mite accession, provided by the URL variable

    Returns:
        The mite entry using the entry.html page as string.
    """

    def count_entries() -> int:
        dirpath = Path(__file__).parent.parent.joinpath("data/data")
        return sum(1 for entry in os.scandir(dirpath) if entry.is_file())

    src = current_app.config["DATA_HTML"].joinpath(f"{mite_acc}.json")

    if not src.exists():
        return render_template("entry_not_found.html", mite_acc=mite_acc)

    with open(src) as infile:
        data = json.load(infile)

    acc_int = int(mite_acc.split("MITE")[1])
    ret_acc = f"MITE{str(acc_int - 1).zfill(7)}"
    fwd_acc = f"MITE{str(acc_int + 1).zfill(7)}"

    max_entry_reached = True if count_entries() == acc_int else False

    return render_template(
        "entry.html",
        data=data,
        ret_acc=ret_acc,
        fwd_acc=fwd_acc,
        max_entry_reached=max_entry_reached,
    )
