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
from datetime import date
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
from pydantic import BaseModel, Field, ValidationError

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
        """

        nr_auxenzymes = set()
        reactions = {}

        for key in data:
            if match := re.search(pattern=r"auxenzyme\[(\d+)\]", string=key):
                nr_auxenzymes.add(match.group(1))
            elif match := re.search(
                pattern=r"reaction\[(\d+)\]knownreaction\[(\d+)\]substrate", string=key
            ):
                if match.group(1) in reactions:
                    reactions[match.group(1)].add(match.group(2))
                else:
                    reactions[match.group(1)] = {match.group(2)}

        nr_auxenzymes = sorted(nr_auxenzymes)

        reactions = {
            "tailoring": ["Macrolactam formation"],
            "description": "Lassopeptide macrolactam formation after N-terminal cleavage: NH2-G-G-X-X-X-X-X-E-Xn-COOH.",
            "reactionSMARTS": "[#6:55]-[#6:54](-[#7:56]-[#6:57](=[O:58])-[#6:59]-[#7:60]-[#6:61](=[O:62])-[#6:63]-[#7:64])-[#6:52](=[O:53])-[#7:51]-[#6:50]-[#6:48](=[O:49])-[#7:47]-[#6@@H:40](-[#6:41]-[c:42]1[c:43][n:44][c:45][n:46]1)-[#6:38](=[O:39])-[#7:37]-[#6@@H:36](-[#6:65])-[#6:34](=[O:35])-[#7:33]-[#6:32]-[#6:30](=[O:31])-[#7:29]-[#6@@H:23](-[#6:24]-[#6:25]-[#6:26](-[#8:28])=[O:27])-[#6:21](=[O:22])-[#7:20]-[#6@@H:18](-[#6:19])-[#6:16](=[O:17])-[#7:15]-[#6@@H:13](-[#6:14])-[#6:11](=[O:12])-[#7:10]-[#6:9]-[#6:7](=[O:8])-[#7:6]-[#6:5]-[#6:3](=[O:4])-[#7:2]-[#6:1]-[#6:66](=[O:67])-[#7:68]-[#6:69]-[#6:70](=[O:71])-[#7:72]-[#6:73]-[#6:74](=[O:75])-[#7:76]-[#6@@H:77](-[#6:111])-[#6:78](=[O:79])-[#7:80]-[#6@@H:81](-[#6:110])-[#6:82](=[O:83])-[#7:84]-[#6@@H:85](-[#6:86])-[#6:87](=[O:88])-[#7:89]-[#6@@H:90](-[#6:91])-[#6:92](=[O:93])-[#7:94]-[#6@@H:95](-[#6:96]-[c:97]1[c:98][c:99][c:100][c:101][c:102]1)-[#6:103](=[O:104])-[#7:105]-[#6:106]-[#6:107](-[#8:109])=[O:108]>>[#6:111]-[#6@H:77](-[#7:76]-[#6:74](=[O:75])-[#6:73]-[#7:72]-[#6:70](=[O:71])-[#6:69]-[#7:68]-[#6:66](=[O:67])-[#6:1]-[#7:2]-[#6:3](=[O:4])-[#6:5]-[#7:6]-[#6:7](=[O:8])-[#6:9]-[#7:10]-[#6:11](=[O:12])-[#6@H:13](-[#6:14])-[#7:15]-[#6:16](=[O:17])-[#6@H:18](-[#6:19])-[#7:20]-[#6:21](=[O:22])-[#6@@H:23]-1-[#6:24]-[#6:25]-[#6:26](=[O:27])-[#7:64]-[#6:63]-[#6:61](=[O:62])-[#7:60]-[#6:59]-[#6:57](=[O:58])-[#7:56]-[#6:54](-[#6:55])-[#6:52](=[O:53])-[#7:51]-[#6:50]-[#6:48](=[O:49])-[#7:47]-[#6@@H:40](-[#6:41]-[c:42]2[c:43][n:44][c:45][n:46]2)-[#6:38](=[O:39])-[#7:37]-[#6@@H:36](-[#6:65])-[#6:34](=[O:35])-[#7:33]-[#6:32]-[#6:30](=[O:31])-[#7:29]-1)-[#6:78](=[O:79])-[#7:80]-[#6@@H:81](-[#6:110])-[#6:82](=[O:83])-[#7:84]-[#6@@H:85](-[#6:86])-[#6:87](=[O:88])-[#7:89]-[#6@@H:90](-[#6:91])-[#6:92](=[O:93])-[#7:94]-[#6@@H:95](-[#6:96]-[c:97]1[c:98][c:99][c:100][c:101][c:102]1)-[#6:103](=[O:104])-[#7:105]-[#6:106]-[#6:107](-[#8:109])=[O:108]",
            "databaseIds": {"rhea": "12345", "ec": "1.2.3.4"},
            "reactions": [
                {
                    "substrate": "CC[C@H](C)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)CNC(=O)[C@H](C)NC(=O)CNC(=O)CN)C(C)C)C(C)C)C(=O)NCC(=O)N[C@H](C(=O)N1CCC[C@H]1C(=O)N[C@H](C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)O)[C@@H](C)CC)[C@@H](C)O",
                    "products": [
                        "CC[C@H](C)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@@H]1CCC(=O)NCC(=O)NCC(=O)N[C@@H](C)C(=O)NCC(=O)N[C@@H](Cc2c[nH]cn2)C(=O)N[C@@H](C(C)C)C(=O)N2CCC[C@H]2C(=O)N1)C(C)C)C(=O)NCC(=O)N[C@H](C(=O)N1CCC[C@H]1C(=O)N[C@H](C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)O)[C@@H](C)CC)[C@@H](C)O"
                    ],
                    "forbidden_products": ["asdfas"],
                    "isIntermediate": False,
                    "description": "Lasso macrocyclisation of microcin J25 after precursor cleavage, leading to mature product.",
                }
            ],
            "evidence": {
                "evidenceCode": [
                    "Site-directed mutagenesis",
                    "Heterologous expression",
                ],
                "references": ["doi:10.1074/jbc.M803995200"],
            },
        }

        self.data = {
            "accession": data["mite_accession"][0]
            if data["mite_accession"][0] != ""
            else "MITE9999999",
            "status": "active",
            "changelog": [],
            "enzyme": {
                "name": data["enzyme_name"][0],
                "description": data["enzyme_description"][0],
                "references": list(data["enzyme_ref[]"]),
                "databaseIds": {
                    "uniprot": data["enzyme_uniprot"][0],
                    "genpept": data["enzyme_uniprot"][0],
                    "mibig": data["enzyme_uniprot"][0],
                },
            },
            "reactions": [],
        }

        for index in nr_auxenzymes:
            if not self.data["enzyme"].get("auxiliaryEnzymes"):
                self.data["enzyme"]["auxiliaryEnzymes"] = []
            self.data["enzyme"]["auxiliaryEnzymes"].append(
                {
                    "name": data.get(f"auxenzyme[{index}]name")[0],
                    "description": data.get(f"auxenzyme[{index}]description")[0],
                    "databaseIds": {
                        "uniprot": data.get(f"auxenzyme[{index}]uniprot")[0],
                        "genpept": data.get(f"auxenzyme[{index}]uniprot")[0],
                    },
                }
            )

        if original_data != {}:
            for version in original_data["changelog"]:
                self.data["changelog"].append(version)

        self.data["changelog"].append(
            {
                "version": f"{len(self.data["changelog"]) + 1}.0",
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

        print(self.data)

    def validate_user_input(self: Self):
        pass


#     raises all the errors
#     catch errors with a try: except in the route, use flash to alert user, return template with the user-data so that it is not necessary to redo the whole page


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
        processing_helper.parse_user_input(data=user_input, original_data=data)

        return user_input

        try:
            processing_helper.validate_user_input()
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
        processing_helper.parse_user_input(data=user_input, original_data={})

        try:
            processing_helper.validate_user_input()
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
