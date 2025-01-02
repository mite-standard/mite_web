"""Generate data files for entry html pages

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
import logging
import shutil
from pathlib import Path

from mite_extras import MiteParser
from pydantic import BaseModel

logger = logging.getLogger("prep_data")


class HtmlJsonManager(BaseModel):
    """Prepare mite jsons for use in html templates

    Attributes:
        src: the location in which the downloaded mite json files are stored
        target: the location in which the converted mite json files will be stored
    """

    src: Path = Path(__file__).parent.parent.joinpath("data/data")
    target: Path = Path(__file__).parent.parent.joinpath("data/data_html")

    def run(self) -> None:
        """Call methods for preparation of html-compatible json files"""
        logger.info("HtmlJsonManager: Started")

        if self.target.exists():
            logger.warning(
                f"DownloadManager: Directory {self.target} already exists - skip json preparation"
            )
            return

        self.target.mkdir(parents=True, exist_ok=True)
        self.convert_json_html()

        logger.info("HtmlJsonManager: Completed")

    def convert_json_html(self) -> None:
        """Convert regular mite json files to html-compatible ones"""
        for entry in self.src.iterdir():
            with open(entry) as infile:
                mite_data = json.load(infile)

            parser = MiteParser()
            parser.parse_mite_json(data=mite_data)

            with open(self.target.joinpath(entry.name), "w") as outfile:
                outfile.write(
                    json.dumps(parser.to_html(), indent=2, ensure_ascii=False)
                )
