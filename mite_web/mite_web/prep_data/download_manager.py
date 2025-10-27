"""Download data from Zenodo and organize

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
import os
import shutil
from pathlib import Path

import requests
from pydantic import BaseModel

logger = logging.getLogger("prep_data")


class DownloadManager(BaseModel):
    """Download data and prepare for use by mite_web

    Attributes:
        mite_data: Zenodo URL for mite_data: always resolves to latest version
        mite_web_extras: Zenodo URL for mite_web_extras: always resolves to latest version
        location: the location to download data to
    """

    mite_data: str = "https://zenodo.org/api/records/13294303"
    mite_web_extras: str = "https://zenodo.org/api/records/17453501"
    location: Path = Path(__file__).parent.parent.joinpath("data")

    def run(self) -> None:
        """Call methods for downloading and moving data"""

        logger.info("DownloadManager: Started")

        if self.location.exists():
            logger.warning(
                f"DownloadManager: Download location {self.location} already exists - skip download"
            )
            return
        self.location.mkdir()

        self.download_mite_data()
        self.organize_mite_data()

        self.download_mite_web_extras()
        self.organize_mite_web_extras()

        self.generate_summary()

        logger.info("DownloadManager: Completed")

    def download_mite_data(self) -> None:
        """Download and store mite_data record from Zenodo

        Raises:
            RuntimeError: Could not download files
        """
        logger.info("DownloadManager: Start mite_data download")

        rsps_meta = requests.get(self.mite_data)
        if rsps_meta.status_code != 200:
            raise RuntimeError(
                f"Error fetching 'mite_data' metadata: {rsps_meta.status_code}"
            )
        record_metadata = rsps_meta.json()
        version = record_metadata["metadata"]["version"]
        files_url = record_metadata["files"][0]["links"]["self"]

        rsps_data = requests.get(files_url)
        if rsps_data.status_code != 200:
            raise RuntimeError(
                f"Error downloading 'mite_data' record: {rsps_data.status_code}"
            )

        with open(self.location.joinpath("mite_data.zip"), "wb") as f:
            f.write(rsps_data.content)

        with open(self.location.joinpath("version.json"), "w") as f:
            f.write(json.dumps({"version_mite_data": f"{version}"}))

        logger.info("DownloadManager: Completed mite_data download")

    def download_mite_web_extras(self) -> None:
        """Download and store mite_web_extras record from Zenodo

        Raises:
            RuntimeError: Could not download files
        """
        logger.info("DownloadManager: Started mite_web_extras download")

        rsps_meta = requests.get(self.mite_web_extras)
        if rsps_meta.status_code != 200:
            raise RuntimeError(
                f"Error fetching 'mite_web_extras' metadata: {rsps_meta.status_code}"
            )
        record_metadata = rsps_meta.json()
        files_url = record_metadata["files"][0]["links"]["self"]

        rsps_data = requests.get(files_url)
        if rsps_data.status_code != 200:
            raise RuntimeError(
                f"Error downloading 'mite_web_extras' record: {rsps_data.status_code}"
            )

        with open(self.location.joinpath("mite_web_extras.zip"), "wb") as f:
            f.write(rsps_data.content)

        logger.info("DownloadManager: Completed mite_web_extras download")

    def organize_mite_data(self) -> None:
        """Unpacks mite_data, moves to correct locations, cleans up

        Raises:
            NotADirectoryError: directory not unzipped in expected location
            RuntimeError: Could not determine data location in downloaded folder
        """
        logger.info("DownloadManager: Started to organize mite_data download")

        src = self.location.joinpath("mite_data.zip")
        trgt = self.location.joinpath("mite_data_unzip")
        download = self.location.joinpath("download")

        download.mkdir(exist_ok=True)

        shutil.unpack_archive(filename=src, extract_dir=trgt, format="zip")
        if not trgt.exists():
            raise NotADirectoryError(f"Could not find the unzipped directory {trgt}.")

        matching_dirs = list(trgt.glob("mite-standard-mite_data-*"))
        if not matching_dirs:
            raise RuntimeError(
                f"Could not determine data storage location in downloaded directory."
            )
        subdir = matching_dirs[0]

        shutil.move(
            src=trgt.joinpath(subdir).joinpath("mite_data/data").resolve(),
            dst=self.location.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir).joinpath("mite_data/fasta").resolve(),
            dst=self.location.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir)
            .joinpath("mite_data/metadata/dump_smarts.csv")
            .resolve(),
            dst=download.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir)
            .joinpath("mite_data/metadata/dump_smiles.csv")
            .resolve(),
            dst=download.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir)
            .joinpath("mite_data/metadata/metadata_general.json")
            .resolve(),
            dst=self.location.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir)
            .joinpath("mite_data/metadata/metadata_mibig.json")
            .resolve(),
            dst=self.location.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir)
            .joinpath("mite_data/metadata/product_list.pickle")
            .resolve(),
            dst=self.location.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir)
            .joinpath("mite_data/metadata/reaction_fps.pickle")
            .resolve(),
            dst=self.location.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir)
            .joinpath("mite_data/metadata/substrate_list.pickle")
            .resolve(),
            dst=self.location.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir)
            .joinpath("mite_data/metadata/summary.csv")
            .resolve(),
            dst=self.location.resolve(),
        )

        os.remove(src)
        shutil.rmtree(trgt)

        logger.info("DownloadManager: Completed to organize mite_data download")

    def organize_mite_web_extras(self) -> None:
        """Unpacks mite_web_extras, moves to correct locations, cleans up

        Raises:
            NotADirectoryError: directory not unzipped in expected location
            RuntimeError: Could not determine data location in downloaded folder
        """
        logger.info("DownloadManager: Started to organize mite_web_extras download")

        src = self.location.joinpath("mite_web_extras.zip")
        trgt = self.location.joinpath("mite_web_extras_unzip")
        download = self.location.joinpath("download")

        download.mkdir(exist_ok=True)

        shutil.unpack_archive(filename=src, extract_dir=trgt, format="zip")
        if not trgt.exists():
            raise NotADirectoryError(f"Could not find the unzipped directory {trgt}.")

        matching_dirs = list(trgt.glob("mite-standard-mite_web_extras-*"))
        if not matching_dirs:
            raise RuntimeError(
                f"Could not determine data storage location in downloaded directory."
            )
        subdir = matching_dirs[0]

        shutil.move(
            src=trgt.joinpath(subdir).joinpath("data/blastlib").resolve(),
            dst=self.location.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir).joinpath("data/html").resolve(),
            dst=self.location.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir).joinpath("data/img").resolve(),
            dst=self.location.parent.joinpath("static").resolve(),
        )
        shutil.move(
            src=trgt.joinpath(subdir).joinpath("data/mite_concat.fasta").resolve(),
            dst=download.resolve(),
        )

        with open(self.location.joinpath("version.json")) as infile:
            v_mite_data = json.load(infile)

        with open(trgt.joinpath(subdir).joinpath("data/version.json")) as infile:
            v_mite_web_extr = json.load(infile)

        if v_mite_data["version_mite_data"] != v_mite_web_extr["version_mite_data"]:
            raise RuntimeError(
                f"Version mite_data {v_mite_data["version_mite_data"]} does not match version mite_web_extras {v_mite_web_extr["version_mite_data"]}. \n"
                f"Was mite_web_extras updated?"
            )

        os.remove(src)
        shutil.rmtree(trgt)

        logger.info(
            "DownloadManager: Completed organization of mite_web_extras download"
        )

    def generate_summary(self) -> None:
        """Creates summary file from metadata_general.json"""
        summary = {"entries": {}}

        with open(self.location.joinpath("metadata_general.json")) as infile:
            metadata = json.load(infile)

        for key, val in metadata["entries"].items():
            summary["entries"][key] = {
                "accession": val["accession"],
                "status": val["status_icon"],
                "status_plain": val["status"],
                "name": val["enzyme_name"],
                "tailoring": val["tailoring"],
                "cofactors_organic": val["cofactors_organic"],
                "cofactors_inorganic": val["cofactors_inorganic"],
                "description": val["enzyme_description"],
                "reaction_description": val["reaction_description"],
                "organism": val["organism"],
                "domain": val["domain"],
                "kingdom": val["kingdom"],
                "phylum": val["phylum"],
                "class": val["class"],
                "order": val["order"],
                "family": val["family"],
            }

        keys = sorted(summary.get("entries").keys())
        summary = {"entries": {key: summary["entries"][key] for key in keys}}

        with open(self.location.joinpath("summary.json"), "w") as f:
            f.write(json.dumps(summary))
