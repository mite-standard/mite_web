"""Routes for main and auxiliary pages.

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
from pathlib import Path

from flask import Response, abort, current_app, jsonify, render_template, send_file

from mite_web.routes import bp


@bp.route("/")
def index() -> str:
    """Render the index page of mite_web

    Returns:
        The index.html page as string.
    """
    return render_template("index.html")


@bp.route("/about/")
def about() -> str:
    """Render the about page of mite_web

    Returns:
        The repository.html page as string.
    """
    return render_template("about.html")


@bp.route("/contact/")
def contact() -> str:
    """Render the contact page of mite_web

    Returns:
        The contact.html page as string.
    """
    return render_template("contact.html")


@bp.route("/terms/")
def termsofuse() -> str:
    """Render the terms of use page of mite_web

    Returns:
        The submission_terms_of_use.html page as string.
    """
    return render_template("submission_terms_of_use.html")


@bp.route("/help/")
def help() -> str:
    """Render the help page of mite_web

    Returns:
        The help.html page as string.
    """
    return render_template("help.html")


@bp.route("/downloads/")
def downloads() -> str:
    """Render the downloads page of mite_web

    Returns:
        The downloads.html page as string.
    """
    return render_template("downloads.html")


@bp.route("/api/<mite_acc>", methods=["GET"])
def api(mite_acc: str) -> Response:
    """Return the MITE entry JSON upon an API call

    Arguments:
        mite_acc: A valid MITE accession (regexp: ^MITE\d\d\d\d\d\d\d$)

    Returns:
        A Response object (a MITE JSON or a 404)
    """
    location = Path(__file__).parent.parent.joinpath(f"data/data/{mite_acc}.json")

    if not location.exists():
        return abort(404, description="Item not found")
    else:
        with open(location) as infile:
            mite = json.load(infile)
            return jsonify(mite)


@bp.route("/downloads/<identifier>")
def download_identifier(identifier: str) -> Response | None:
    """Delivers the file in question for download to user

    Arguments:
        identifier: the string identifier of the file

    Returns:
        A Response object containing the file or None
    """
    if identifier == "smarts":
        return send_file(
            Path(__file__).parent.parent.joinpath("data/download/dump_smarts.csv"),
            as_attachment=True,
        )
    elif identifier == "smiles":
        return send_file(
            Path(__file__).parent.parent.joinpath("data/download/dump_smiles.csv"),
            as_attachment=True,
        )
    elif identifier == "blastlib":
        return send_file(
            Path(__file__).parent.parent.joinpath("data/download/MiteBlastDB.zip"),
            as_attachment=True,
        )
    elif identifier == "mite_zip":
        return send_file(
            Path(__file__).parent.parent.joinpath(
                "data/download/MITE_all_active_entries.zip"
            ),
            as_attachment=True,
        )
    elif identifier.startswith("MITE"):
        return send_file(
            Path(__file__).parent.parent.joinpath(f"data/data/{identifier}.json"),
            as_attachment=True,
        )
    elif identifier == "mite_fasta":
        return send_file(
            Path(__file__).parent.parent.joinpath("data/download/mite_concat.fasta"),
            as_attachment=True,
        )
    else:
        return
