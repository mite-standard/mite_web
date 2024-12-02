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
import pickle
from pathlib import Path

import pandas as pd
from flask import current_app, flash, render_template, request
from rdkit.Chem import MolFromSmarts, MolFromSmiles, PandasTools

from mite_web.routes import bp


@bp.route("/overview/", methods=["GET", "POST"])
def overview() -> str:
    """Render the repository overview page of mite_web

    Returns:
        The overview.html page as string.
    """

    with open(current_app.config["DATA_SUMMARY"]) as infile:
        summary = json.load(infile)

    if request.method == "POST":
        user_input = request.form.to_dict()

        if (
            not user_input.get("substructure_query")
            or user_input.get("substructure_query") == ""
        ):
            flash("Please specify a substructure query string.")
            return render_template("overview.html", entries=summary.get("entries"))

        try:
            df = pd.read_csv(
                Path(__file__).parent.parent.joinpath("data/download/dump_smiles.csv")
            )

            with open(
                Path(__file__).parent.parent.joinpath("data/substrate_list.pickle"),
                "rb",
            ) as infile:
                substrate_list = pickle.load(infile)

            with open(
                Path(__file__).parent.parent.joinpath("data/product_list.pickle"), "rb"
            ) as infile:
                product_list = pickle.load(infile)

            df["ROMol_substrates"] = substrate_list
            df["ROMol_products"] = product_list

            if user_input.get("action") == "smiles":
                mol_query = MolFromSmiles(user_input.get("substructure_query"))
            else:
                mol_query = MolFromSmarts(user_input.get("substructure_query"))

            df_substrate_match = df[df["ROMol_substrates"] >= mol_query]
            df_product_match = df[df["ROMol_products"] >= mol_query]

            unique_mite_acc = set(df_substrate_match["mite_id"].str.split(".").str[0])
            unique_mite_acc.update(
                set(df_product_match["mite_id"].str.split(".").str[0])
            )

            matching_entries = {
                key: value
                for key, value in summary.get("entries").items()
                if key in unique_mite_acc
            }

            return render_template(
                "overview.html",
                entries=matching_entries,
                query_data={
                    "substructure_query": user_input.get("substructure_query"),
                    "result_len": len(matching_entries),
                },
            )

        except Exception as e:
            flash(f"An error in the substructure matching occurred: '{e!s}'")
            return render_template("overview.html", entries=summary.get("entries"))

    return render_template("overview.html", entries=summary.get("entries"))


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
