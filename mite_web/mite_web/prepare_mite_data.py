"""Populates the app with data from mite_data

This script runs as part of the Dockerfile during container-building.
It automatically downloads the latest version of mite_data from Zenodo,
generates the html-json files, and prepares metadata.

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
import shutil
import sys
from pathlib import Path

import requests
from mite_extras import MiteParser
from pydantic import BaseModel


class DownloadManager(BaseModel):
    """Download data and prepare for use by mite_web

    Attributes:
        record_url: Zenodo URL for mite_data: always resolves to latest version
        location: the location to download data to
        record: path to the record file
        record_unzip: path to unzipped record file
        version: path to the mite_data version location
    """

    record_url: str = "https://zenodo.org/api/records/13294303"
    location: Path = Path(__file__).parent.joinpath("data")
    record: Path = Path(__file__).parent.joinpath("data/record.zip")
    record_unzip: Path = Path(__file__).parent.joinpath("data/record")
    version: Path = Path(__file__).parent.joinpath("data/version.json")

    def run(self) -> None:
        """Call methods for downloading and moving data"""
        if self.location.exists():
            print(
                f"Location '{self.location}' exists already - data is presumably already downloaded. SKIP"
            )
            return

        self.location.mkdir(parents=True)
        self.download_data()
        self.organize_data()

    def download_data(self) -> None:
        """Download data from Zenodo

        Raises:
            RuntimeError: Could not download files
        """
        response_metadata = requests.get(self.record_url)
        if response_metadata.status_code != 200:
            raise RuntimeError(
                f"Error fetching 'mite_data' record metadata: {response_metadata.status_code}"
            )

        record_metadata = response_metadata.json()
        version = record_metadata["metadata"]["version"]
        files_url = record_metadata["files"][0]["links"]["self"]

        response_data = requests.get(files_url)
        if response_data.status_code != 200:
            raise RuntimeError(
                f"Error downloading 'mite_data' record: {response_data.status_code}"
            )

        with open(self.record, "wb") as f:
            f.write(response_data.content)
        with open(self.version, "w") as f:
            f.write(json.dumps({"version_mite_data": f"{version}"}))

    def organize_data(self) -> None:
        """Unpacks data, moves to convenient location, cleans up

        Raises:
            NotADirectoryError: directory not unzipped in expected location
        """
        shutil.unpack_archive(
            filename=self.record, extract_dir=self.record_unzip, format="zip"
        )

        if not self.record_unzip.exists():
            raise NotADirectoryError(
                f"Could not find the unzipped directory {self.record_unzip}."
            )

        matching_dirs = list(self.record_unzip.glob("mite-standard-mite_data-*"))
        if not matching_dirs:
            raise RuntimeError(
                f"Could not determine data storage location in downloaded directory."
            )

        subdir = matching_dirs[0]

        shutil.move(
            src=self.record_unzip.joinpath(subdir).joinpath("mite_data/data").resolve(),
            dst=self.location.resolve(),
        )
        shutil.move(
            src=self.record_unzip.joinpath(subdir).joinpath("mite_data/img").resolve(),
            dst=self.location.parent.joinpath("static/").resolve(),
        )
        shutil.move(
            src=self.record_unzip.joinpath(subdir)
            .joinpath("mite_data/blast_lib")
            .resolve(),
            dst=self.location.resolve(),
        )

        os.remove(self.record)
        shutil.rmtree(self.record_unzip)


class HtmlJsonManager(BaseModel):
    """Prepare mite jsons for use in html templates

    Attributes:
        src: the location in which the downloaded mite json files are stored
        target: the location in which the converted mite json files will be stored
    """

    src: Path = Path(__file__).parent.joinpath("data/data")
    target: Path = Path(__file__).parent.joinpath("data/data_html")

    def run(self) -> None:
        """Call methods for preparation of html-compatible json files"""
        if self.target.exists():
            print(
                f"Location '{self.target}' exists already - data is presumably already converted. SKIP"
            )
            return

        self.target.mkdir(parents=True)
        self.convert_json_html()

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


class AuxFileManager(BaseModel):
    """Prepare auxiliary files for website

    Attributes:
        src: the location in which the downloaded mite json files are stored
        target: the location in which the auxiliary files will be stored
    """

    src: Path = Path(__file__).parent.joinpath("data/data")
    target: Path = Path(__file__).parent.joinpath("data/")

    def run(self) -> None:
        """Call methods for preparation of auxiliary files"""
        self.prepare_summary()

    def prepare_summary(self) -> None:
        """Create a summary of mite entries for repository table"""
        summary = {"entries": {}}

        for entry in self.src.iterdir():
            with open(entry) as infile:
                mite_data = json.load(infile)

            reviewer = set()
            for log in mite_data.get("changelog"):
                reviewer.update(log.get("reviewers"))

            summary["entries"][mite_data.get("accession")] = {
                "status": '<i class="bi bi-check-circle-fill"></i>'
                if mite_data.get("status") == "active"
                else '<i class="bi bi-circle"></i>',
                "reviewed": '<i class="bi bi-check-circle-fill"></i>'
                if reviewer != {"BBBBBBBBBBBBBBBBBBBBBBBB"}
                else '<i class="bi bi-circle"></i>',
                "name": mite_data.get("enzyme", {}).get("name"),
                "tailoring": "|".join(
                    sorted(
                        {
                            tailoring
                            for reaction in mite_data.get("reactions")
                            for tailoring in reaction.get("tailoring", [])
                        }
                    )
                ),
                "description": mite_data.get("enzyme", {}).get("description", "N/A"),
                "reaction_description": mite_data["reactions"][0].get(
                    "description", "No description available"
                ),
            }

        keys = list(summary.get("entries").keys())
        keys.sort()
        summary_sorted = {"entries": {key: summary["entries"][key] for key in keys}}

        with open(self.target.joinpath("summary.json"), "w") as outfile:
            outfile.write(json.dumps(summary_sorted, indent=2, ensure_ascii=False))


def main() -> None | SystemExit:
    """Prepares mite_data for use in mite_web"""
    try:
        download_manager = DownloadManager()
        download_manager.run()

        json_manager = HtmlJsonManager()
        json_manager.run()

        aux_manager = AuxFileManager()
        aux_manager.run()
    except Exception as e:
        print(f"An exception has occurred: {e}")
        return sys.exit(1)


if __name__ == "__main__":
    main()
