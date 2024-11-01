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

from flask import Response, current_app, redirect, render_template, request, url_for
from mite_extras.processing.data_classes import EnyzmeDatabaseIds
from mite_extras.processing.mite_parser import MiteParser
from pydantic import BaseModel, Field, ValidationError

from mite_web.routes import bp


@bp.route("/submission/")
def submission() -> str:
    """Render the submission page of mite_web

    Returns:
        The submission.html page as string.
    """
    return render_template("submission.html")


@bp.route("/submission/<mite_acc>/", methods=["GET", "POST"])
def submission_existing(mite_acc: str) -> str | Response:
    """Render the submission forms for an existing entry mite_acc

    Arguments:
        mite_acc: the mite accession, provided by the URL variable

    Returns:
        The submission_existing.html page as string.
    """
    if request.method == "POST":
        user_input = request.form.to_dict(flat=False)
        return user_input

    src = current_app.config["DATA_JSON"].joinpath(f"{mite_acc}.json")

    if not src.exists():
        return render_template("entry_not_found.html", mite_acc=mite_acc)

    with open(src) as infile:
        data = json.load(infile)

    if data["status"] != "active":
        return redirect(url_for("routes.repository", mite_acc=mite_acc))

    print(data)

    # strategy: use the json directly for form creation
    # makro: take name etc + value; can be Null to not set a default value; this way, reusable
    # not all data is needed; better to specify manually
    # we are not going to majorly change the schema

    return render_template("submission_form.html", data=data)
