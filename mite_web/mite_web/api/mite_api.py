"""Adds rudimentary API for MITE

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
import re
from pathlib import Path

from flask import abort, jsonify
from flask_restx import Namespace, Resource

mite_ns = Namespace("mite", description="Operations related to MITE files")


@mite_ns.route("/<string:mite_acc>")
@mite_ns.param(
    "mite_acc",
    "A MITE accession in the format: 'MITE' followed by seven digits (e.g., 'MITE0000001')",
)
@mite_ns.response(404, "Item not found.")
@mite_ns.response(
    400, "Invalid MITE accession format. Must be 'MITE' followed by seven digits."
)
@mite_ns.response(200, "Successful retrieval of MITE JSON entry.")
class MITEApi(Resource):
    def get(self, mite_acc):
        """Retrieve a MITE JSON entry.
        Returns a MITE JSON file or an appropriate error message.
        """
        if not re.fullmatch(r"MITE\d{7}", mite_acc):
            abort(
                400,
                description="Invalid MITE accession format. Must be 'MITE' followed by seven digits.",
            )

        location = Path(__file__).parent.parent.joinpath(f"data/data/{mite_acc}.json")
        if not location.exists():
            abort(404, description="Item not found")

        with open(location) as infile:
            mite = json.load(infile)
            return jsonify(mite)
