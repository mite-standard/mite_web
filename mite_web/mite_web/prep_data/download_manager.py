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
        record_url: Zenodo URL for mite_data: always resolves to latest version
        location: the location to download data to
        record: path to the record file
        record_unzip: path to unzipped record file
        version: path to the mite_data version location
    """

    record_url: str = "https://zenodo.org/api/records/13294303"
    location: Path = Path(__file__).parent.parent.joinpath("data")
    record: Path = Path(__file__).parent.parent.joinpath("data/record.zip")
    record_unzip: Path = Path(__file__).parent.parent.joinpath("data/record")
    version: Path = Path(__file__).parent.parent.joinpath("data/version.json")

    def run(self) -> None:
        """Call methods for downloading and moving data"""

        logger.info("DownloadManager: Started")

        if self.location.exists():
            logger.warning(
                f"DownloadManager: Download location {self.location} already exists - skip download"
            )
            return

        self.location.mkdir(parents=True)
        self.download_data()
        self.organize_data()

        logger.info("DownloadManager: Completed")

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
            RuntimeError: Could not determine data location in downloaded folder
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
            src=self.record_unzip.joinpath(subdir)
            .joinpath("mite_data/fasta")
            .resolve(),
            dst=self.location.resolve(),
        )

        # TODO(MMZ 21.12.24): replace with de-novo generation of MITE BLASTLIE
        self.location.joinpath("download/").mkdir(exist_ok=True)
        shutil.copy(
            src=str(
                self.record_unzip.joinpath(subdir).joinpath(
                    "mite_data/blast_lib/MiteBlastDB.zip"
                )
            ),
            dst=str(self.location.joinpath("download/MiteBlastDB.zip")),
        )
        self.location.joinpath("blastlib/").mkdir(exist_ok=True)
        shutil.copy(
            src=str(
                self.record_unzip.joinpath(subdir).joinpath(
                    "mite_data/blast_lib/MiteBlastDB.zip"
                )
            ),
            dst=str(self.location.joinpath("blastlib/MiteBlastDB.zip")),
        )
        shutil.unpack_archive(
            filename=self.location.joinpath("blastlib/MiteBlastDB.zip"),
            extract_dir=self.location.joinpath("blastlib/"),
            format="zip",
        )
        os.remove(self.location.joinpath("blastlib/MiteBlastDB.zip"))

        os.remove(self.record)
        shutil.rmtree(self.record_unzip)
