"""Manages downloading mibig gene references.

Copyright (c) 2024 to present Mitja Maximilian Zdouc, PhD and individual contributors.

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

import logging
import os
from pathlib import Path

import pandas as pd
import requests
from pydantic import BaseModel

logger = logging.getLogger("prep_data")


class MibigManager(BaseModel):
    """Download data and prepare for use by mite_web

    Attributes:
        record_url: Zenodo URL for mibig: always resolves to latest version
        location: the location to download data to
        record: path to the record file
    """

    record_url: str = "https://zenodo.org/api/records/13367755"
    location: Path = Path(__file__).parent.parent.joinpath("data")
    record: Path = Path(__file__).parent.parent.joinpath("data/mibig_proteins.fasta")

    def run(self) -> None:
        """Call methods for downloading and moving data"""

        logger.info("MibigManager: Started")

        self.location.mkdir(exist_ok=True)
        self.download_data()
        self.organize_data()

        logger.info("MibigManager: Completed")

    def download_data(self) -> None:
        """Download data from Zenodo

        Raises:
            RuntimeError: Could not download files
        """
        response_metadata = requests.get(self.record_url)
        if response_metadata.status_code != 200:
            raise RuntimeError(
                f"Error fetching 'mibig' record metadata: {response_metadata.status_code}"
            )

        record_metadata = response_metadata.json()

        for entry in record_metadata["files"]:
            if entry["key"].endswith("fasta"):
                response_data = requests.get(entry["links"]["self"])

                if response_data.status_code != 200:
                    raise RuntimeError(
                        f"Error downloading 'mibig' record: {response_data.status_code}"
                    )

                with open(self.record, "wb") as f:
                    f.write(response_data.content)

        if not self.record.exists():
            raise RuntimeError(
                f"Could not find the mibig fasta file in its Zenodo repository (record URL: {self.record_url})"
            )

    def organize_data(self) -> None:
        """Extract data, move to location

        Raises:
            NotADirectoryError: directory not unzipped in expected location
            RuntimeError: Could not determine data location in downloaded folder
        """
        mibig_prot = {
            "mibig": [],
            "genpept": [],
        }

        with open(self.record) as infile:
            for line in infile.readlines():
                if line.startswith(">"):
                    accs = line.split("|")
                    mibig_prot["mibig"].append(accs[0].removeprefix(">").split(".")[0])
                    mibig_prot["genpept"].append(accs[-1].replace("\n", ""))

        df = pd.DataFrame(mibig_prot)
        df.to_csv(self.location.joinpath("mibig_proteins.csv"), index=False)
        os.remove(self.record)
