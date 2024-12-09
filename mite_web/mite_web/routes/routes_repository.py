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

import json
import os
import pickle
import re
import subprocess
import uuid
from pathlib import Path

import pandas as pd
from Bio import Blast, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from flask import current_app, flash, render_template, request
from pydantic import BaseModel
from rdkit.Chem import MolFromSmarts, MolFromSmiles, PandasTools

from mite_web.routes import bp


class QueryManager(BaseModel):
    """Organize querying functions

    Attributes:
        summary: a dict of the entries
        dump_smiles: path to csv-file of all SMILES in MITE entries
        substrates_pickle: path to pre-calculated fingerprints of substrate SMILES
        products_pickle: path to pre-calculated fingerprints of products SMILES
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
    blastlib: Path = Path(__file__).parent.parent.joinpath("data/blastlib/")

    def return_summary(self) -> dict:
        """Returns the summary dict"""
        return self.summary

    def query_substructure(self, action: str, query: str) -> dict:
        """Query dataset for a specific substructure and filter summary for matching entries

        Arguments:
            action: the type of matching to perform (smiles or smarts matching)
            query: a SMILES or SMARTS string

        Returns:
            A dict with substructure matching results
        """

        df = pd.read_csv(self.dump_smiles)

        with open(
            self.substrates_pickle,
            "rb",
        ) as infile:
            substrate_list = pickle.load(infile)

        with open(self.products_pickle, "rb") as infile:
            product_list = pickle.load(infile)

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

        mite_entries = []
        substructure_results = {}
        for idx, hit in enumerate(matches, 1):
            mite_entries.append(hit.split(".")[0])
            substructure_results[idx] = {
                "accession": hit.split(".")[0],
                "reaction": hit.split(".")[1],
                "example": hit.split(".")[2],
                "query": query,
            }

        self.summary = {
            key: value for key, value in self.summary.items() if key in mite_entries
        }

        return substructure_results

    def query_sequence(self, query: str, e_val: int) -> dict:
        """Run BLAST against MITE BLAST DB for a specific protein and filter summary for matching entries

        Arguments:
            query: a protein sequence
            e_val: the e-value to filter matches with

        Returns:
            A dictionary of BLAST results

        Raises:
            RuntimeError: input sanitization detected illegal input
        """
        query = query.replace("\n", "")
        query = re.sub(r"\s+", "", query, flags=re.UNICODE)

        if not query.isalpha():
            print(query)
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

        blast_results = {}
        mite_entries = []
        for idx, hit in enumerate(blast_record, 1):
            mite_entries.append(hit[0].target.description.split()[0])
            blast_results[idx] = {
                "accession": hit[0].target.description.split()[0],
                "sequence_similarity": round(
                    (
                        (float(hit[0].annotations.get("positive")) / float(len(query)))
                        * 100
                    ),
                    0,
                ),
                "alignment_score": round(hit[0].score, 0),
                "evalue": round(hit[0].annotations.get("evalue"), 0),
                "bit_score": round(hit[0].annotations.get("bit score"), 0),
            }
        os.remove(self.blastlib.joinpath(f"{job_uuid}.xml"))
        os.remove(self.blastlib.joinpath(f"{job_uuid}.fasta"))

        self.summary = {
            key: value for key, value in self.summary.items() if key in mite_entries
        }

        return blast_results


@bp.route("/overview/", methods=["GET", "POST"])
def overview() -> str:
    """Render the repository overview page of mite_web

    Returns:
        The overview.html page as string.
    """

    with open(current_app.config["DATA_SUMMARY"]) as infile:
        summary_json = json.load(infile)
        summary = summary_json.get("entries")

    if request.method == "POST":
        query_manager = QueryManager(summary=summary)
        user_input = request.form.to_dict()

        if user_input.get("action") == "smiles" or user_input.get("action") == "smarts":
            if user_input.get("substructure_query") == "":
                flash("Please specify a substructure query string.")
                return render_template("overview.html", entries=summary)

            try:
                results = query_manager.query_substructure(
                    action=user_input.get("action"),
                    query=user_input.get("substructure_query"),
                )
                return render_template(
                    "overview.html",
                    entries=query_manager.return_summary(),
                    substructure_results=results,
                )
            except Exception as e:
                flash(f"An error in the substructure matching occurred: '{e!s}'")
                return render_template("overview.html", entries=summary)

        elif user_input.get("action") == "blast":
            if user_input.get("sequence") == "":
                flash("Please specify a amino acid sequence query string.")
                return render_template("overview.html", entries=summary)

            try:
                blast_results = query_manager.query_sequence(
                    query=user_input.get("sequence"), e_val=int(user_input.get("e_val"))
                )
                return render_template(
                    "overview.html",
                    entries=query_manager.return_summary(),
                    blast_results=blast_results,
                )

            except Exception as e:
                flash(f"An error in BLASTp matching occurred: '{e!s}'")
                return render_template("overview.html", entries=summary)

    return render_template("overview.html", entries=summary)


@bp.route("/repository/<mite_acc>/")
def repository(mite_acc: str) -> str:
    """Render the individual pages

    Arguments:
        mite_acc: the mite accession, provided by the URL variable

    Returns:
        The mite entry using the entry.html page as string.
    """
    src = current_app.config["DATA_HTML"].joinpath(f"{mite_acc}.json")

    if not src.exists():
        return render_template("entry_not_found.html", mite_acc=mite_acc)

    with open(src) as infile:
        data = json.load(infile)

    acc_int = int(mite_acc.split("MITE")[1])
    ret_acc = f"MITE{str(acc_int - 1).zfill(7)}"
    fwd_acc = f"MITE{str(acc_int + 1).zfill(7)}"

    return render_template("entry.html", data=data, ret_acc=ret_acc, fwd_acc=fwd_acc)
