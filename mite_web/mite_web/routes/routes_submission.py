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
import re
import shutil
import subprocess
import time
import uuid
from datetime import date
from pathlib import Path
from typing import Self

import pandas as pd
import requests
from flask import (
    Response,
    abort,
    current_app,
    flash,
    redirect,
    render_template,
    request,
    send_file,
    session,
    url_for,
)
from mite_extras.processing.mite_parser import MiteParser
from mite_schema import SchemaManager
from pydantic import BaseModel
from rdkit import Chem

from mite_web.routes import bp


def get_schema_vals() -> dict:
    """Extract values from MITE JSON Schema"""
    with open(SchemaManager().entry) as f:
        schema = json.load(f)

    return {
        "evidence": schema["$defs"]["evidence"]["enum"],
        "tailoring": schema["$defs"]["tailoringFunction"]["enum"],
        "inorganic": schema["$defs"]["inorganic"]["enum"],
        "organic": schema["$defs"]["organic"]["enum"],
    }


def create_validated_parser(data: dict) -> MiteParser:
    """Creates validated MiteParser instance"""
    parser = MiteParser()
    parser.parse_mite_json(data=data)
    schema_manager = SchemaManager()
    schema_manager.validate_mite(instance=parser.to_json())
    return parser


def render_preview(data: dict) -> dict:
    """Run validations and render html-json"""
    parser = create_validated_parser(data)
    return parser.to_html()


class ProcessingHelper(BaseModel):
    """Contains methods to help processing the data

    Attributes:
        dump_name: a name under which the file is dumped
        data: the user-submitted data
        reviewer_tags: registered reviewers
    """

    dump_name: str
    data: dict | None = None
    # TODO(MMZ 24.7.25): add reviewers after briefing
    reviewer_tags: tuple = ("@mmzdouc",)
    reviewer_orcids: tuple = ("0000-0001-6534-6609",)

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
                    "wikidata": data.get("enzyme_wikidata", [""])[0]
                    if data.get("enzyme_wikidata", [""])[0] != ""
                    else None,
                },
                "cofactors": {
                    "organic": data.get(f"enzyme-cofactors-organic-check[]", []),
                    "inorganic": data.get(f"enzyme-cofactors-inorganic-check[]", []),
                },
            },
            "reactions": [],
        }

        if data.get("comment", [""])[0] != "":
            self.data["comment"] = data.get("comment")[0]

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
                        "wikidata": data.get(f"auxenzyme[{index}]wikidata", [""])[0],
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

    def add_changelog(self, form: dict) -> None:
        """Parses changelog information and adds to existing mite entry changelog

        Arguments:
            form: the user-input
        """
        self.data["changelog"].append(
            {
                "version": f"{len(self.data["changelog"]) + 1}",
                "date": date.today().strftime("%Y-%m-%d"),
                "contributors": [
                    form["orcid"] if form["orcid"] != "" else "AAAAAAAAAAAAAAAAAAAAAAAA"
                ],
                "reviewers": ["BBBBBBBBBBBBBBBBBBBBBBBB"],
                "comment": form["changelog"],
            }
        )

    def add_reviewer_info(self, form: dict) -> None:
        """Parses reviewer information and adds to existing mite entry changelog

        Arguments:
            form: the user-input

        Raises:
            RuntimeError: reviewer orcid not part of allowed orcids
        """
        reviewer = form["reviewer-orcid"]
        if reviewer not in self.reviewer_orcids:
            raise RuntimeError(
                f"ORCID '{reviewer}' is not one of registered reviewer ORCIDs. If you want to become a reviewer for MITE, contact the developers."
            )

        current_reviewers = set(self.data["changelog"][-1]["reviewers"])
        current_reviewers.add(reviewer)
        current_reviewers.discard("BBBBBBBBBBBBBBBBBBBBBBBB")
        self.data["changelog"][-1]["reviewers"] = list(current_reviewers)

        if form["reviewer-changelog"] != "":
            self.data["changelog"][-1]["comment"] = (
                self.data["changelog"][-1]["comment"] + " " + form["reviewer-changelog"]
            )

    def validate_user_input(self: Self) -> None:
        """Validates the incoming user-submitted and formatted data

        Raises:
            RuntimeError: input validation does not pass
        """
        if not self.data["reactions"]:
            raise RuntimeError("Please provide at least one reaction entry!")

        enzyme_db_ids = [
            self.data["enzyme"]["databaseIds"].get("uniprot"),
            self.data["enzyme"]["databaseIds"].get("genpept"),
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

    # TODO(MMZ 23.7): re-implement as webhook
    def check_mibig(self: Self) -> None:
        """Check if genpept ID can be found in mibig genes

        Raises:
            RuntimeError: genpept found in MIBiG protein accession list
        """
        genpept = self.data["enzyme"]["databaseIds"].get("genpept")

        df = pd.read_csv(
            Path(__file__).parent.parent.joinpath("data/mibig_proteins.csv")
        )
        matches = df[df["genpept"].str.contains(genpept)]

        if len(matches) > 0:
            raise RuntimeError(
                f"NCBI GenPept Accession '{genpept}' is associated to MIBiG entry '{matches["mibig"].iloc[0]}', but was not added in this form. Please consider adding this cross-reference. This message will appear only once."
            )

    # TODO(MMZ 23.7): re-implement as webhook
    def check_rhea(self: Self) -> None:
        """Check if rhea ids can be found for uniprot and/or are already added

        Raises:
            RuntimeError: found non-overlapping ids between uniprot rheas and form
        """
        rhea = set()
        form = set()

        response = requests.get(
            url="https://www.rhea-db.org/rhea?",
            params={
                "query": self.data["enzyme"]["databaseIds"].get("uniprot"),
                "columns": "rhea-id",
                "format": "tsv",
                "limit": 10,
            },
            timeout=3,
        )
        if response.status_code == 200:
            rhea = {i.removeprefix("RHEA:") for i in response.text.split()[2:]}

        for reaction in self.data["reactions"]:
            if val := reaction.get("databaseIds", {}).get("rhea"):
                form.add(val)

        if len(rhea) == 0:
            return
        else:
            diff = form.intersection(rhea)
            if len(diff) == 0:
                raise RuntimeError(
                    f"The UniProt ID '{self.data["enzyme"]["databaseIds"].get("uniprot")}' is described in Rhea, but its Rhea IDs were not detected in the submission form. Please consider adding this cross-reference. This message will appear only once. The detected Rhea IDs are: {rhea}."
                )

    def dump_json(self: Self):
        """Dumps dict as JSON to disk"""
        with open(
            current_app.config["DATA_DUMPS"].joinpath(self.dump_name),
            "w",
            encoding="utf-8",
        ) as outfile:
            outfile.write(json.dumps(self.data, indent=4, ensure_ascii=False))

    def create_pr(self) -> None:
        """Create PR on mite_data using mite_bot's credentials"""

        if not current_app.config.get("ONLINE", False):
            current_app.logger.warning(
                f"{self.dump_name}: Prevented PR in offline mode"
            )
            return

        src = current_app.config["DATA_DUMPS"].joinpath(f"{self.dump_name}")
        trgt = current_app.config["MITE_DATA"].joinpath(
            f"mite_data/data/{self.dump_name}"
        )
        if self.data["accession"] != "MITE9999999":
            trgt = current_app.config["MITE_DATA"].joinpath(
                f"mite_data/data/{self.data["accession"]}.json"
            )
        shutil.copy(src, trgt)

        branch = self.dump_name.split(".")[0]

        body = f"""
A submission was performed via the MITE web portal and needs reviewing.

## Review requested

{", ".join(self.reviewer_tags)}

## TODO Reviewers

- Review the entry [HERE](https://mite.bioinformatics.nl/submission/preview/{branch}/reviewer)
- Fix any issues, add your ORCID, download the file, and append it to this PR

*This action was performed by `mite-bot`*
"""
        subprocess.run(
            ["git", "-C", current_app.config["MITE_DATA"], "checkout", "-b", branch],
            check=True,
        )
        subprocess.run(
            ["git", "-C", current_app.config["MITE_DATA"], "add", "mite_data/data/"],
            check=True,
        )
        subprocess.run(
            [
                "git",
                "-C",
                current_app.config["MITE_DATA"],
                "commit",
                "-m",
                f"Contributor submission {self.data['enzyme']['name']}",
            ],
            check=True,
        )

        subprocess.run(
            [
                "git",
                "-C",
                current_app.config["MITE_DATA"],
                "push",
                "-u",
                "origin",
                branch,
            ],
            check=True,
        )

        subprocess.run(
            [
                "gh",
                "pr",
                "create",
                "--repo",
                "mite-standard/mite_data",
                "--title",
                f"Contributor submission {branch}",
                "--body",
                "Automated submission from webapp",
                "--draft",
                "--base",
                "main",
                "--head",
                branch,
                "--body",
                body,
            ],
            check=True,
        )

        subprocess.run(
            ["git", "-C", current_app.config["MITE_DATA"], "checkout", "main"],
            check=True,
        )
        subprocess.run(
            ["git", "-C", current_app.config["MITE_DATA"], "pull"], check=True
        )
        subprocess.run(
            ["git", "-C", current_app.config["MITE_DATA"], "branch", "-D", branch],
            check=True,
        )


@bp.route("/submission/")
def submission() -> str:
    """Render the submission page of mite_web

    Returns:
        The submission.html page as string.
    """
    return render_template("submission.html")


@bp.route("/submission/<var>/<role>", methods=["GET", "POST"])
def submission_data(var: str, role: str) -> str | Response:
    """Initiate new submission and render form page

    Arguments:
        var: 'new', a mite accession id, or an uuid -> determine which type of data
        role: 'contributor' or 'reviewer' to indicate the role of data modification

    Returns:
        Form page or redirect to 'entry_not_found' or 'retired' pages
    """
    if var.startswith("MITE"):
        src = current_app.config["DATA_JSON"].joinpath(f"{var}.json")
        if not src.exists():
            return render_template("entry_not_found.html", mite_acc=var)

        with open(src) as infile:
            data = json.load(infile)
            if data.get("status") == "retired":
                return redirect(url_for("routes.repository", mite_acc=var))

        var = uuid.uuid1()
        with open(
            current_app.config["DATA_DUMPS"].joinpath(f"{var}.json"),
            "w",
            encoding="utf-8",
        ) as h:
            h.write(json.dumps(data, indent=4, ensure_ascii=False))

    elif var == "new":
        data = {
            "changelog": [],
            "enzyme": {"references": ["doi:"]},
            "reactions": [
                {
                    "evidence": {"evidenceCode": [], "references": ["doi:"]},
                    "reactions": [{"products": [""]}],
                }
            ],
        }
        var = uuid.uuid1()
        with open(
            current_app.config["DATA_DUMPS"].joinpath(f"{var}.json"),
            "w",
            encoding="utf-8",
        ) as h:
            h.write(json.dumps(data, indent=4, ensure_ascii=False))

    else:
        if not current_app.config["DATA_DUMPS"].joinpath(f"{var}.json").exists():
            return redirect(url_for("routes.submission_data", var="new"))

    with open(current_app.config["DATA_DUMPS"].joinpath(f"{var}.json")) as infile:
        data = json.load(infile)

    return render_template(
        "submission_form.html",
        data=data,
        form_vals=get_schema_vals(),
        var=var,
        role=role,
    )


@bp.route("/submission/process/<var>/<role>", methods=["GET", "POST"])
def submission_process(var: str, role: str) -> str | Response:
    """Process submitted data for preview

    Arguments:
        var: an uuid to determine storage location
        role: 'contributor' or 'reviewer' to indicate the role of data modification

    Returns:
        If submission error, show rendered page; else redirect to preview page
    """
    if not current_app.config["DATA_DUMPS"].joinpath(f"{var}.json").exists():
        return redirect(
            url_for("routes.submission_data", var="new", role="contributor")
        )

    if request.method == "POST":
        user_input = request.form.to_dict(flat=False)
        processing_helper = ProcessingHelper(dump_name=f"{var}.json")

        with open(current_app.config["DATA_DUMPS"].joinpath(f"{var}.json")) as infile:
            data = json.load(infile)

        try:
            processing_helper.parse_user_input(data=user_input, original_data=data)
        except Exception as e:
            return render_template("submission_failure.html", error=str(e))

        try:
            processing_helper.validate_user_input()
            processing_helper.dump_json()
            return redirect(url_for("routes.submission_preview", var=var, role=role))
        except Exception as e:
            processing_helper.dump_json()
            current_app.logger.critical(
                f"{var}: Error during validation of submission: {e!s}"
            )
            flash(str(e))
            return redirect(url_for("routes.submission_data", var=var, role=role))

    else:
        return redirect(
            url_for("routes.submission_data", var="new", role="contributor")
        )


@bp.route("/submission/preview/<var>/<role>", methods=["GET", "POST"])
def submission_preview(var: str, role: str) -> str | Response:
    """Render the preview page

    preview=True triggers preview class in base.html, adding preview watermark

    Arguments:
        var: the submission ID
        role: the client role (contributor or reviewer)

    Returns:
        The entry.html page as string or entry not found or submission_success
    """
    src = current_app.config["DATA_DUMPS"].joinpath(f"{var}.json")

    if not src.exists():
        return render_template("entry_not_found.html", mite_acc=var)

    with open(src) as infile:
        data = json.load(infile)

    if request.method == "POST":
        user_input = request.form.to_dict()
        processing_helper = ProcessingHelper(dump_name=f"{var}.json", data=data)

        if user_input.get("contr-submit"):
            if request.form.get("email_confirm"):
                # This is a bot (field should be empty)
                abort(400)

            elapsed = time.time() - session.get("form_start", 0)
            if elapsed < 2:  # less than 2 seconds? likely a bot
                abort(400)

            processing_helper.add_changelog(user_input)
            processing_helper.dump_json()
            processing_helper.create_pr()

            return render_template(
                "submission_success.html", sub_id=Path(processing_helper.dump_name).stem
            )
        elif user_input.get("contr-modify"):
            return redirect(
                url_for("routes.submission_data", var=var, role="contributor")
            )
        elif user_input.get("reviewer-modify"):
            return redirect(url_for("routes.submission_data", var=var, role="reviewer"))
        elif user_input.get("reviewer-submit"):
            try:
                processing_helper.add_reviewer_info(user_input)
            except RuntimeError as e:
                flash(str(e))
                current_app.logger.error(f"{e!s}")
                return render_template(
                    "entry.html",
                    data=render_preview(data),
                    mode="review",
                    preview=True,
                    submission_id=var,
                )
            processing_helper.dump_json()
            return send_file(
                current_app.config["DATA_DUMPS"].joinpath(f"{var}.json"),
                as_attachment=True,
            )

        # todo: condition reviewer download -> returns the download as dump with added orcid or reviewer

    if role == "contributor":
        session["form_start"] = time.time()
        return render_template(
            "entry.html",
            data=render_preview(data),
            mode="preview",
            preview=True,
            submission_id=var,
        )
    else:
        return render_template(
            "entry.html",
            data=render_preview(data),
            mode="review",
            preview=True,
            submission_id=var,
        )


# TODO: change the URL to only peptidesmiles
@bp.route("/submission/peptidesmiles", methods=["GET", "POST"])
def peptidesmiles() -> str:
    """Render the peptide SMILES page

    Returns:
        The smiles_peptide.html page as string.
    """
    if request.method == "POST":
        user_input = request.form.to_dict()
        try:
            return render_template(
                "smiles_peptide.html",
                data={
                    "peptide_string": user_input.get("peptide_string"),
                    "smiles": f"{Chem.MolToSmiles(Chem.MolFromSequence(user_input.get("peptide_string")))}",
                },
            )
        except Exception as e:
            current_app.logger.critical(e)
            flash(str(e))
            return render_template(
                "smiles_peptide.html",
                data={"peptide_string": user_input.get("peptide_string", "")},
            )

    return render_template("smiles_peptide.html", data={"peptide_string": ""})


# TODO: change the URL to only canonicalizesmiles
@bp.route("/submission/canonicalizesmiles", methods=["GET", "POST"])
def canonsmiles() -> str:
    """Render the canonicalize SMILES page

    Returns:
        The smiles_canonicalize.html page as string.
    """
    if request.method == "POST":
        user_input = request.form.to_dict()
        try:
            return render_template(
                "smiles_canonicalize.html",
                data={
                    "smiles_in": user_input.get("smiles_in"),
                    "smiles_out": f"{Chem.MolToSmiles(Chem.MolFromSmiles(user_input.get("smiles_in")))}",
                },
            )
        except Exception as e:
            current_app.logger.critical(e)
            flash(str(e))
            return render_template(
                "smiles_canonicalize.html",
                data={"smiles_in": user_input.get("smiles_in", "")},
            )

    return render_template("smiles_canonicalize.html", data={"smiles_in": ""})
