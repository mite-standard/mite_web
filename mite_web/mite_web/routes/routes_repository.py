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
from typing import Any
import uuid
from pathlib import Path

import pandas as pd
from Bio import Blast, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from flask import current_app, flash, render_template, request
from pydantic import BaseModel
from rdkit.Chem import MolFromSmarts, MolFromSmiles, PandasTools, rdChemReactions
from rdkit.DataStructs import FingerprintSimilarity
from sqlalchemy import or_, and_
from sqlalchemy.orm import class_mapper

from mite_web.routes import bp
from mite_web.models import Entry, Enzyme, ChangeLog, Cofactor

class QueryManager(BaseModel):
    """Organize querying functions

    Attributes:
        summary: a dict of the entries
        dump_smiles: path to csv-file of all SMILES in MITE entries
        substrates_pickle: path to pre-calculated fingerprints of substrate SMILES
        products_pickle: path to pre-calculated fingerprints of products SMILES
        reaction_pickle: path to pre-calculated fingerprints of reaction SMARTS
        blastlib: path to the blast library
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
    blastlib: Path = Path(__file__).parent.parent.joinpath("data/blastlib/")

    def return_summary(self) -> dict:
        """Returns the summary dict

        Returns:
            Modified summary dict

        Raises:
            RuntimeError: No matches found
        """
        if len(self.summary) == 0:
            raise RuntimeError("No significant matches found")
        else:
            return self.summary

    def query_substructure(self, action: str, query: str):
        """Query dataset for a specific substructure and filter summary for matching entries

        Arguments:
            action: the type of matching to perform (smiles or smarts matching)
            query: a SMILES or SMARTS string
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

        for match in matches:
            key = match.split(".")[0]
            reaction = match.split(".")[1].replace("reaction", "rx")
            example = match.split(".")[2].replace("example", "ex")

            if self.summary[key].get("reaction"):
                self.summary[key]["reaction"].append(reaction)
                self.summary[key]["example"].append(f"{reaction}-{example}")
            else:
                self.summary[key]["reaction"] = [reaction]
                self.summary[key]["example"] = [f"{reaction}-{example}"]

        copy_summary = copy.deepcopy(self.summary)
        for key, value in copy_summary.items():
            if value["status"] == '<i class="bi bi-circle"></i>' or not value.get(
                "reaction"
            ):
                self.summary.pop(key, None)
            else:
                self.summary[key]["reaction"] = ", ".join(
                    [str(i) for i in self.summary[key]["reaction"]]
                )
                self.summary[key]["example"] = ", ".join(
                    [str(i) for i in self.summary[key]["example"]]
                )

    def query_reaction(self, query: str, sim_score: float):
        """Query dataset for a reaction and filter summary for matching entries

        Arguments:
            query: a reaction SMARTS string
            sim_score: the minimum similarity score to consider an entry
        """

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

        df = pd.DataFrame(reaction_dict)
        for _, row in df.iterrows():
            key = row["mite_id"].split(".")[0]
            if row["similarities"] >= sim_score:
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

        copy_summary = copy.deepcopy(self.summary)
        for key, value in copy_summary.items():
            if value["status"] == '<i class="bi bi-circle"></i>' or not value.get(
                "reaction"
            ):
                self.summary.pop(key, None)
            else:
                self.summary[key]["reaction"] = ", ".join(
                    [str(i) for i in self.summary[key]["reaction"]]
                )
                self.summary[key]["sim_score"] = ", ".join(
                    [str(i) for i in self.summary[key]["sim_score"]]
                )

    def query_sequence(self, query: str, e_val: int):
        """Run BLAST against MITE BLAST DB for a specific protein and filter summary for matching entries

        Arguments:
            query: a protein sequence
            e_val: the e-value to filter matches with

        Raises:
            RuntimeError: input sanitization detected illegal input
        """
        query = query.replace("\n", "")
        query = re.sub(r"\s+", "", query, flags=re.UNICODE)

        if not query.isalpha():
            raise RuntimeError("The query AA sequence contains illegal characters.")
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

        for hit in blast_record:
            key = hit[0].target.description.split()[0]
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

        copy_summary = copy.deepcopy(self.summary)
        for key, value in copy_summary.items():
            if value["status"] == '<i class="bi bi-circle"></i>' or not value.get(
                "sequence_similarity"
            ):
                self.summary.pop(key, None)



class DatabaseManager:
    """Organize functions for database querying

    Attributes:
        operators: translates querybuilder operators to SQLAlchemy col, val pairs
    """
    operators : dict = {
        'equal': lambda col, val: col == val,
        'not_equal': lambda col, val: col != val,
        'contains': lambda col, val: col.ilike(f"%{val}%"),
        'not_contains': lambda col, val: ~col.ilike(f"%{val}%"),
        'is_null': lambda col, _: col.is_(None),
        'is_not_null': lambda col, _: col.is_not(None),
    }

    field_map : dict = {
        "accession": "Entry.accession",
        "contributor": "Entry.changelogs.contributors.orcid",
        "reviewer": "Entry.changelogs.reviewers.orcid",


        "enzyme.mibig_id": "Entry.enzyme.mibig_id",
        "enzyme.uniprot_id": "Entry.enzyme.uniprot_id",
        "enzyme.cofactors.cofactor_name": "Entry.enzyme.cofactors.cofactor_name"
    }

    def query_db(self, rules: dict) -> set:
        """Queries the DB

        Args:
            rules: the raw query from QueryBuilder JSON rules

        Returns:
            A set of MITE accession IDs for filtering
        """
        query = Entry.query
        filters, join_paths = self.parse_rules(rules)

        for path in join_paths:
            current_query = query
            for rel_attr in path:
                current_query = current_query.join(rel_attr)
            query = current_query

        if filters is not None:
            query = query.filter(filters)

        return {e.accession for e in query.all()}

    def parse_rules(self, rules: dict) -> Any | None:
        """Converts QueryBuilder JSON rules into a SQLAlchemy-compatible expression

        Args:
            rules: QueryBuilder JSON rules

        Returns:
            A SQLAlchemy filter expression
        """
        filters = []
        join_paths = set()

        for rule in rules["rules"]:
            field_path = self.field_map.get(rule["field"])
            if not field_path:
                continue

            column, joins = self.resolve_field_path(field_path)
            if joins:
                join_paths.add(tuple(joins))

            operator = rule["operator"]
            value = rule.get('value')

            if column is not None and operator in self.operators:
                filters.append(self.operators[operator](column, value))

        if not filters:
            return None, join_paths

        condition = rules.get("condition", "AND").upper()
        combined = and_(*filters) if condition == "AND" else or_(*filters)
        return combined, join_paths


    def resolve_field_path(self, field_path: str) -> tuple[Any, list]:
        """Convert string path into column, join_path

        Args:
            field_path: a mapping from a QueryBuilder rule

        Returns:
            A tuple of a column name and a list of relationship attributes to use in query.join()
        """
        parts = field_path.split(".")
        current_model = globals()[parts[0]]
        join_path = []
        attr = None

        for part in parts[1:]:
            mapper = class_mapper(current_model)

            if part in mapper.relationships:
                join_path.append(getattr(current_model, part))
                current_model = mapper.relationships[part].mapper.class_
            else:
                attr = getattr(current_model, part)

        return attr, join_path

@bp.route("/overview/", methods=["GET", "POST"])
def overview() -> str:
    """Render the repository overview page of mite_web

    Returns:
        The overview.html page as string.
    """
    summary = current_app.config["SUMMARY"]
    accessions = current_app.config["ACCESSIONS"]


    if request.method == "POST":

        forms = request.form.to_dict()

        rules = json.loads(request.form['rules'])

        if len(rules["rules"]) > 0:
            db_manager = DatabaseManager()
            accessions.intersection_update(db_manager.query_db(rules))

        # TODO: implement rest of filters in the same way as with dbmanager




    # if request.method == "POST":
    #     query_manager = QueryManager(summary=summary)
    #     user_input = request.form.to_dict()
    #
    #     if user_input.get("action") == "smiles" or user_input.get("action") == "smarts":
    #         if user_input.get("substructure_query") == "":
    #             flash("Please specify a substructure query string.")
    #             return render_template("overview.html", entries=summary)
    #
    #         try:
    #             query_manager.query_substructure(
    #                 action=user_input.get("action"),
    #                 query=user_input.get("substructure_query"),
    #             )
    #             return render_template(
    #                 "overview.html",
    #                 entries=query_manager.return_summary(),
    #             )
    #         except Exception as e:
    #             flash(f"An error in the substructure matching occurred: '{e!s}'")
    #             return render_template("overview.html", entries=summary)
    #
    #     if user_input.get("action") == "reaction":
    #         if user_input.get("reaction_query") == "":
    #             flash("Please specify a reaction SMARTS query.")
    #             return render_template("overview.html", entries=summary)
    #
    #         try:
    #             query_manager.query_reaction(
    #                 query=user_input.get("reaction_query"),
    #                 sim_score=float(user_input.get("similarity")),
    #             )
    #             return render_template(
    #                 "overview.html",
    #                 entries=query_manager.return_summary(),
    #             )
    #         except Exception as e:
    #             flash(f"An error in the reaction matching occurred: '{e!s}'")
    #             return render_template("overview.html", entries=summary)
    #
    #     elif user_input.get("action") == "blast":
    #         if user_input.get("sequence") == "":
    #             flash("Please specify a amino acid sequence query string.")
    #             return render_template("overview.html", entries=summary)
    #
    #         try:
    #             query_manager.query_sequence(
    #                 query=user_input.get("sequence"), e_val=int(user_input.get("e_val"))
    #             )
    #             return render_template(
    #                 "overview.html",
    #                 entries=query_manager.return_summary(),
    #             )
    #         except Exception as e:
    #             flash(f"An error in BLASTp matching occurred: '{e!s}'")
    #             return render_template("overview.html", entries=summary)

    summary = {key: val for key, val in summary.items() if key in accessions}

    return render_template("overview.html", entries=summary)


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
