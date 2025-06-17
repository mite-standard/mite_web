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
from io import BytesIO
from pathlib import Path
from typing import Self

import pandas as pd
import requests
from flask import (
    Response,
    current_app,
    flash,
    redirect,
    render_template,
    request,
    url_for,
)
from flask_mail import Message
from mite_extras.processing.mite_parser import MiteParser
from mite_schema import SchemaManager
from pydantic import BaseModel
from rdkit import Chem

from mite_web.config.extensions import mail
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


class ProcessingHelper(BaseModel):
    """Contains methods to help processing the data

    Attributes:
        dump_name: a name under which the file is dumped
        data: the user-submitted data
    """

    dump_name: str
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
                    "wikidata": data.get("enzyme_wikidata", [""])[0]
                    if data.get("enzyme_wikidata", [""])[0] != ""
                    else None,
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

        self.data["changelog"].append(
            {
                "version": f"{len(self.data["changelog"]) + 1}",
                "date": date.today().strftime("%Y-%m-%d"),
                "contributors": [
                    data["orcid"][0]
                    if data["orcid"][0] != ""
                    else "AAAAAAAAAAAAAAAAAAAAAAAA"
                ],
                "reviewers": ["BBBBBBBBBBBBBBBBBBBBBBBB"],
                "comment": data["changelog"][0],
            }
        )

    def validate_user_input(self: Self, initial: str):
        """Validates the incoming user-submitted and formatted data

        Arguments:
            initial: a string ("true", "false") indicating if validation has failed previously

        Raises:
            RuntimeError: input validation does not pass
        """
        if not self.data["reactions"]:
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
            if re.search(r"\|", reaction.get("reactionSMARTS")):
                raise RuntimeError(
                    f"Reaction SMARTS with CXSMARTS (Chemaxon SMARTS) elements detected which are not supported by MITE. The offending reaction SMARTS is: '{reaction.get("reactionSMARTS")}'"
                )

        if initial == "true":
            if (
                self.data["enzyme"]["databaseIds"].get("genpept")
                and not self.data["enzyme"]["databaseIds"]["mibig"]
            ):
                self.check_mibig()

            if self.data["enzyme"]["databaseIds"].get("uniprot") and not self.data[
                "enzyme"
            ]["databaseIds"].get("uniprot").startswith("UPI"):
                self.check_rhea()

        parser = MiteParser()
        parser.parse_mite_json(data=self.data)

        schema_manager = SchemaManager()
        schema_manager.validate_mite(instance=parser.to_json())

        self.data = parser.to_json()

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
        target = Path(__file__).parent.parent.joinpath("dumps")
        target.mkdir(parents=True, exist_ok=True)

        with open(target.joinpath(self.dump_name), "w", encoding="utf-8") as outfile:
            outfile.write(json.dumps(self.data, indent=4, ensure_ascii=False))

    def send_email(self: Self, sub_type: str) -> None:
        """Sends data per email to data processor

        Arguments:
            sub_type: a string indicating the type of submission (modification of entry or new entry)
        """

        if current_app.config.get("ONLINE", False):
            msg = Message()
            msg.recipients = [current_app.config.get("MAIL_TARGET")]
            msg.subject = f"MITE: {sub_type}"
            msg.body = f"Please find the dump file '{self.dump_name}' attached."

            json_content = json.dumps(self.data, indent=4)
            json_attachment = BytesIO(json_content.encode("utf-8"))
            msg.attach(self.dump_name, "application/json", json_attachment.read())

            try:
                mail.send(msg)
                current_app.logger.info(
                    f"Data of file {self.dump_name} was emailed successfully"
                )
            except Exception as e:
                current_app.logger.error(
                    f"An error occurred during sending of data of file {self.dump_name}: {e}."
                )


@bp.route("/submission/")
def submission() -> str:
    """Render the submission page of mite_web

    Returns:
        The submission.html page as string.
    """
    return render_template("submission.html")


@bp.route("/submission/<mite_acc>", methods=["GET", "POST"])
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

        processing_helper = ProcessingHelper(dump_name=f"{uuid.uuid1()}.json")

        try:
            processing_helper.parse_user_input(data=user_input, original_data=data)
        except Exception as e:
            return render_template("submission_failure.html", error=str(e))

        try:
            processing_helper.validate_user_input(initial=user_input["initial"][0])
            processing_helper.dump_json()
            processing_helper.send_email(sub_type="MODIFIED")
            return render_template(
                "submission_success.html", sub_id=Path(processing_helper.dump_name).stem
            )
        except Exception as e:
            current_app.logger.critical(e)
            flash(str(e))
            x, y = processing_helper.random_numbers()
            return render_template(
                "submission_form.html",
                data=processing_helper.data,
                x=x,
                y=y,
                initial="false",
                form_vals=get_schema_vals(),
            )

    x, y = ProcessingHelper(dump_name=f"{uuid.uuid1()}.json").random_numbers()

    return render_template(
        "submission_form.html",
        data=data,
        x=x,
        y=y,
        initial="true",
        form_vals=get_schema_vals(),
    )


@bp.route("/submission/new", methods=["GET", "POST"])
def submission_new() -> str | Response:
    """Render the submission forms for a new entry

    Returns:
        The submission_existing.html page as string or redirect to another page
    """
    if request.method == "POST":
        user_input = request.form.to_dict(flat=False)

        processing_helper = ProcessingHelper(dump_name=f"{uuid.uuid1()}.json")

        try:
            processing_helper.parse_user_input(data=user_input, original_data={})
        except Exception as e:
            return render_template("submission_failure.html", error=str(e))

        try:
            processing_helper.validate_user_input(initial=user_input["initial"][0])
            processing_helper.dump_json()
            processing_helper.send_email(sub_type="NEW")
            return render_template(
                "submission_success.html", sub_id=Path(processing_helper.dump_name).stem
            )
        except Exception as e:
            current_app.logger.critical(e)
            flash(str(e))
            x, y = processing_helper.random_numbers()
            return render_template(
                "submission_form.html",
                data=processing_helper.data,
                x=x,
                y=y,
                initial="false",
                form_vals=get_schema_vals(),
            )

    x, y = ProcessingHelper(dump_name=f"{uuid.uuid1()}.json").random_numbers()
    data = {
        "enzyme": {"references": ["doi:"]},
        "reactions": [
            {
                "evidence": {"evidenceCode": [], "references": ["doi:"]},
                "reactions": [{"products": [""]}],
            }
        ],
    }

    return render_template(
        "submission_form.html",
        data=data,
        x=x,
        y=y,
        initial="true",
        form_vals=get_schema_vals(),
    )


@bp.route("/submission/review", methods=["GET", "POST"])
def review() -> str:
    """Render the review page

    Returns:
        The entry.html page as string
    """

    def _render(data: dict) -> dict:
        """Run validations and render html-json"""
        parser = MiteParser()
        parser.parse_mite_json(data=data)
        schema_manager = SchemaManager()
        schema_manager.validate_mite(instance=parser.to_json())
        return parser.to_html()

    if request.method == "POST":
        try:
            json_data = ""
            user_input = request.form.to_dict(flat=True)

            if user_input.get("jsonText") and user_input.get("jsonText") != "":
                json_data = json.loads(user_input.get("jsonText"))
            elif "jsonFile" in request.files:
                json_file = request.files["jsonFile"]
                if json_file.filename.endswith(".json"):
                    json_data = json.load(json_file)
            else:
                raise RuntimeError(
                    f"Neither MITE file nor content provided. Please try again."
                )

            return render_template("entry.html", data=_render(json_data))

        except Exception as e:
            current_app.logger.critical(e)
            flash(str(e))
            return render_template("review.html")

    return render_template("review.html")


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
