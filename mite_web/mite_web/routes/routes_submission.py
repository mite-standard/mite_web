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
import random
import re
import uuid
from datetime import date
from pathlib import Path
from typing import Self

from flask import (
    Response,
    current_app,
    flash,
    redirect,
    render_template,
    request,
    url_for,
)
from mite_extras.processing.mite_parser import MiteParser
from mite_schema import SchemaManager
from pydantic import BaseModel

from mite_web.routes import bp


class ProcessingHelper(BaseModel):
    """Contains methods to help processing the data

    Attributes:
        data: the user-submitted data
    """

    data: dict | None = None

    @staticmethod
    def random_numbers() -> tuple[int, int]:
        """Generate two random numbers and return them

        Returns:
            A tuple of two single-digit numbers
        """
        return random.randint(1, 9), random.randint(1, 9)

    def parse_user_input(self: Self, data: dict, original_data: dict):
        """Reads the user_input json dict and brings it in the mite-format

        Arguments:
            data: the user-input
            original_data: the currently established version of the entry/data

        Raises:
            RuntimeError: input validation does not pass
        """
        nr_auxenzymes = set()
        reactions = {}

        for string in data:
            if match := re.search(pattern=r"auxenzyme\[(\d+)\]", string=string):
                nr_auxenzymes.add(match.group(1))
            elif match := re.search(pattern=r"reaction\[(\d+)\]", string=string):
                if match.group(1) not in reactions:
                    reactions[match.group(1)] = set()

        for key in reactions:
            for string in data:
                if match := re.search(
                    pattern=rf"reaction\[({key})\]knownreaction\[(\d+)\]", string=string
                ):
                    reactions[key].add(match.group(2))

        self.data = {
            "accession": data["mite_accession"][0]
            if data["mite_accession"][0] != ""
            else "MITE9999999",
            "status": "pending",
            "changelog": [],
            "enzyme": {
                "name": data.get("enzyme_name", [""])[0],
                "description": data.get("enzyme_description", [""])[0],
                "references": data.get("enzyme_ref[]", []),
                "databaseIds": {
                    "uniprot": data.get("enzyme_uniprot", [""])[0]
                    if data.get("enzyme_uniprot", [""])[0] != ""
                    else None,
                    "genpept": data.get("enzyme_genpept", [""])[0]
                    if data.get("enzyme_genpept", [""])[0] != ""
                    else None,
                    "mibig": data.get("enzyme_mibig", [""])[0]
                    if data.get("enzyme_mibig", [""])[0] != ""
                    else None,
                },
            },
            "reactions": [],
        }

        nr_auxenzymes = sorted(nr_auxenzymes)
        for index in nr_auxenzymes:
            if not self.data["enzyme"].get("auxiliaryEnzymes"):
                self.data["enzyme"]["auxiliaryEnzymes"] = []
            self.data["enzyme"]["auxiliaryEnzymes"].append(
                {
                    "name": data.get(f"auxenzyme[{index}]name", [""])[0],
                    "description": data.get(f"auxenzyme[{index}]description", [""])[0],
                    "databaseIds": {
                        "uniprot": data.get(f"auxenzyme[{index}]uniprot", [""])[0],
                        "genpept": data.get(f"auxenzyme[{index}]genpept", [""])[0],
                    },
                }
            )

        for key, value in reactions.items():
            self.data["reactions"].append(
                {
                    "tailoring": data.get(f"reaction[{key}]tailoring[]", [""]),
                    "description": data.get(f"reaction[{key}]description", [""])[0],
                    "reactionSMARTS": data.get(f"reaction[{key}]smarts", [""])[0],
                    "databaseIds": {
                        "rhea": data.get(f"reaction[{key}]rhea", [""])[0],
                        "ec": data.get(f"reaction[{key}]ec", [""])[0],
                    },
                    "evidence": {
                        "evidenceCode": data.get(
                            f"reaction[{key}]evidencecode[]", [""]
                        ),
                        "references": data.get(f"reaction[{key}]ref[]", [""]),
                    },
                    "reactions": [],
                }
            )
            sorted_val = sorted(value)
            for index in sorted_val:
                self.data["reactions"][int(key)]["reactions"].append(
                    {
                        "substrate": data.get(
                            f"reaction[{key}]knownreaction[{index}]substrate", [""]
                        )[0],
                        "products": data.get(
                            f"reaction[{key}]knownreaction[{index}]products[]", [""]
                        ),
                        "forbidden_products": data.get(
                            f"reaction[{key}]knownreaction[{index}]forbiddenproducts[]"
                        ),
                        "isIntermediate": data.get(
                            f"reaction[{key}]knownreaction[{index}]intermediate", [""]
                        )[0],
                        "description": data.get(
                            f"reaction[{key}]knownreaction[{index}]description", [""]
                        )[0],
                    }
                )

        if original_data != {}:
            for version in original_data["changelog"]:
                self.data["changelog"].append(version)

        self.data["changelog"].append(
            {
                "version": f"{len(self.data["changelog"]) + 1}",
                "date": date.today().strftime("%Y-%m-%d"),
                "contributors": [
                    data["orcid"][0]
                    if data["orcid"][0] != ""
                    else "AAAAAAAAAAAAAAAAAAAAAAAA"
                ],
                "reviewers": ["AAAAAAAAAAAAAAAAAAAAAAAA"],
                "comment": data["changelog"][0],
            }
        )

    def validate_user_input(self: Self):
        """Validates the incoming user-submitted and formatted data

        Raises:
            RuntimeError: input validation does not pass
        """
        if self.data.get("reaction[0]smarts", [""]) == [""]:
            raise RuntimeError("Please provide at least one reaction entry!")

        enzyme_db_ids = [
            self.data["enzyme"]["databaseIds"]["uniprot"],
            self.data["enzyme"]["databaseIds"]["genpept"],
        ]
        if all(item is None for item in enzyme_db_ids):
            raise RuntimeError(
                "No enzyme database cross-references specified. Please provide at least an UniProt or GenPept ID."
            )

        if len(self.data.get("enzyme").get("references")) == 0 or all(
            item == "" for item in self.data.get("enzyme").get("references")
        ):
            raise RuntimeError(
                "No enzyme references specified. Please provide at least one reference. Please contact the Developers if you cannot provide a reference but still want to submit your data to MITE (e.g. because the work is not published yet)."
            )

        for reaction in self.data.get("reactions"):
            if len(reaction.get("evidence", {}).get("references")) == 0 or all(
                item == "" for item in reaction.get("evidence", {}).get("references")
            ):
                raise RuntimeError(
                    "No reaction references specified. Please provide at least one reference. Please contact the Developers if you cannot provide a reference but still want to submit your data to MITE (e.g. because the work is not published yet)."
                )
            if reaction.get("evidence", {}).get("evidenceCode") == [""]:
                raise RuntimeError(
                    "At least one of the checkboxes in 'Experimental Evidence Qualifiers' must be checked."
                )
            if reaction.get("tailoring") == [""]:
                raise RuntimeError(
                    "At least one of the checkboxes in 'Tailoring Reaction Controlled Vocabulary' must be checked."
                )

        parser = MiteParser()
        parser.parse_mite_json(data=self.data)

        schema_manager = SchemaManager()
        schema_manager.validate_mite(instance=parser.to_json())

        self.data = parser.to_json()

    def dump_json(self: Self):
        """Dumps dict as JSON to disk"""
        target = Path(__file__).parent.parent.joinpath("dumps")
        target.mkdir(parents=True, exist_ok=True)

        with open(
            target.joinpath(f"{uuid.uuid1()}.json"), "w", encoding="utf-8"
        ) as outfile:
            outfile.write(json.dumps(self.data, indent=4, ensure_ascii=False))


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
        The submission_existing.html page as string or a redirect to another page
    """
    src = current_app.config["DATA_JSON"].joinpath(f"{mite_acc}.json")

    if not src.exists():
        return render_template("entry_not_found.html", mite_acc=mite_acc)

    with open(src) as infile:
        data = json.load(infile)

    if data.get("status") != "active":
        return redirect(url_for("routes.repository", mite_acc=mite_acc))

    if request.method == "POST":
        user_input = request.form.to_dict(flat=False)

        processing_helper = ProcessingHelper()

        try:
            processing_helper.parse_user_input(data=user_input, original_data=data)
        except Exception as e:
            return render_template("submission_failure.html", error=str(e))

        try:
            processing_helper.validate_user_input()
            processing_helper.dump_json()
            # TODO(MMZ 8.11.24): add emailing step
            return redirect(url_for("routes.submission_success"))
        except Exception as e:
            current_app.logger.critical(e)
            flash(str(e))
            x, y = processing_helper.random_numbers()
            return render_template(
                "submission_form.html", data=processing_helper.data, x=x, y=y
            )

    x, y = ProcessingHelper().random_numbers()

    return render_template("submission_form.html", data=data, x=x, y=y)


@bp.route("/submission/new/", methods=["GET", "POST"])
def submission_new() -> str | Response:
    """Render the submission forms for a new entry

    Returns:
        The submission_existing.html page as string or redirect to another page
    """
    if request.method == "POST":
        user_input = request.form.to_dict(flat=False)

        processing_helper = ProcessingHelper()

        try:
            processing_helper.parse_user_input(data=user_input, original_data={})
        except Exception as e:
            return render_template("submission_failure.html", error=str(e))

        try:
            processing_helper.validate_user_input()
            processing_helper.dump_json()
            # TODO(MMZ 7.11): send data via email and/or dump for testing
            return redirect(url_for("routes.submission_success"))
        except Exception as e:
            current_app.logger.critical(e)
            flash(str(e))
            x, y = processing_helper.random_numbers()
            return render_template(
                "submission_form.html", data=processing_helper.data, x=x, y=y
            )

    x, y = ProcessingHelper().random_numbers()
    data = {
        "enzyme": {"references": ["doi:"]},
        "reactions": [
            {
                "evidence": {"evidenceCode": [], "references": ["doi:"]},
                "reactions": [{"products": [""]}],
            }
        ],
    }

    return render_template("submission_form.html", data=data, x=x, y=y)


@bp.route("/submission/success/")
def submission_success() -> str:
    """Render the successful submission page

    Returns:
        The submission_success.html page as string.
    """
    return render_template("submission_success.html")
