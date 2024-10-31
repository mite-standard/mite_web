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
from pathlib import Path

from flask import current_app, render_template, request
from pydantic import BaseModel, EmailStr, Field, ValidationError

from mite_web.routes import bp


@bp.route("/submission/")
def submission() -> str:
    """Render the submission page of mite_web

    Returns:
        The submission.html page as string.
    """
    return render_template("submission.html")


class UserModel(BaseModel):
    username: str = Field(default="MMM", min_length=3, max_length=50)
    email: EmailStr
    age: int = Field(default=10, gt=0)


@bp.route("/submission/<mite_acc>/", methods=["GET", "POST"])
def submission_existing(mite_acc: str) -> str:
    """Render the submission forms for an existing entry mite_acc

    Arguments:
        mite_acc: the mite accession, provided by the URL variable

    Returns:
        The submission_existing.html page as string.
    """

    # check if the entry exists in the files, if not, leave
    # load the data of the entry
    # put it into the pydantic model

    #####

    # import a subpart, e.g. changelog
    # hard code data
    # try to

    if request.method == "POST":
        print(request.form)
        data = request.form.to_dict()
        print(data)
        try:
            user = UserModel(**data)
            return "Form submitted successfully!"
        except ValidationError as e:
            return str(e), 400

    user_schema = UserModel.model_json_schema()
    return render_template("submission_existing.html", user_schema=user_schema)
