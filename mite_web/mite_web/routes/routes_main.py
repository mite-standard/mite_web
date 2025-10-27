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

import re
import uuid

from flask import Response, current_app, render_template, send_file

from mite_web.routes import bp


@bp.route("/")
def index() -> str:
    """Render the index page of mite_web

    Returns:
        The index.html page as string.
    """
    return render_template("index.html")


@bp.route("/about")
def about() -> str:
    """Render the about page of mite_web

    Returns:
        The repository.html page as string.
    """
    return render_template("about.html")


@bp.route("/contact")
def contact() -> str:
    """Render the contact page of mite_web

    Returns:
        The contact.html page as string.
    """
    return render_template("contact.html")


@bp.route("/terms")
def termsofuse() -> str:
    """Render the terms of use page of mite_web

    Returns:
        The submission_terms_of_use.html page as string.
    """
    return render_template("submission_terms_of_use.html")


@bp.route("/tutorial")
def tutorial() -> str:
    """Render the tutorial page of mite_web

    Returns:
        The tutorial.html page as string.
    """
    return render_template("tutorial.html")


@bp.route("/troubleshooting")
def troubleshooting() -> str:
    """Render the troubleshooting page of mite_web

    Returns:
        The troubleshooting.html page as string.
    """
    return render_template("troubleshooting.html")


@bp.route("/faqs")
def faqs() -> str:
    """Render the faqs page of mite_web

    Returns:
        The faqs.html page as string.
    """
    return render_template("faqs.html")


@bp.route("/downloads/")
def downloads() -> str:
    """Render the downloads page of mite_web

    Returns:
        The downloads.html page as string.
    """
    return render_template("downloads.html")


@bp.route("/downloads/<identifier>")
def download_identifier(identifier: str) -> Response | None:
    """Delivers the file in question for download to user

    'else' condition: assumes a search query uuid ID

    Arguments:
        identifier: the string identifier of the file

    Returns:
        A Response object containing the file or None
    """
    try:
        if identifier == "smarts":
            return send_file(
                current_app.config["DOWNLOAD"].joinpath("dump_smarts.csv"),
                as_attachment=True,
            )
        elif identifier == "smiles":
            return send_file(
                current_app.config["DOWNLOAD"].joinpath("dump_smiles.csv"),
                as_attachment=True,
            )
        elif identifier == "mite_fasta":
            return send_file(
                current_app.config["DOWNLOAD"].joinpath("mite_concat.fasta"),
                as_attachment=True,
            )
        elif identifier == "mite_overview":
            return send_file(
                current_app.config["DOWNLOAD"].parent.joinpath("summary.csv"),
                as_attachment=True,
            )
        elif re.fullmatch(r"MITE[0-9]{7}", identifier):
            path = current_app.config["DATA_JSON"].joinpath(f"{identifier}.json")

            if not path.is_file():
                raise ValueError(f"{identifier} does not exist as MITE entry")
            else:
                return send_file(
                    path,
                    as_attachment=True,
                )
        else:
            uuid.UUID(identifier.split(".")[0])

            if current_app.config["QUERIES"].joinpath(f"{identifier}").is_file():
                return send_file(
                    current_app.config["QUERIES"].joinpath(f"{identifier}"),
                    as_attachment=True,
                )
            elif (
                current_app.config["DATA_DUMPS"]
                .joinpath(f"{identifier}.json")
                .is_file()
            ):
                return send_file(
                    current_app.config["DATA_DUMPS"].joinpath(f"{identifier}.json"),
                    as_attachment=True,
                )
            else:
                raise ValueError(f"{identifier} does not exist")

    except Exception as e:
        current_app.logger.warning(f"Download failed - {e!s}")
        return Response(status=204)
