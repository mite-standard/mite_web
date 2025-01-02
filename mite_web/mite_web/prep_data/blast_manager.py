"""Manages building of a BLAST database.

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

import json
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Self

from pydantic import BaseModel, DirectoryPath

logger = logging.getLogger("prep_data")


class BlastManager(BaseModel):
    """Manages the building of a BLAST database.

    Attributes:
        src: a Path to the souce directory containing fasta files
        target_blast: a Path towards the blast directory (internal use)
        target_download: a Path towards the download directory
        concat_filename: filename of the concatenated fasta file
    """

    src: DirectoryPath = Path(__file__).parent.parent.joinpath("data/fasta/")
    target_download: DirectoryPath = Path(__file__).parent.parent.joinpath(
        "data/download/"
    )
    target_blast: DirectoryPath = Path(__file__).parent.parent.joinpath(
        "data/blastlib/"
    )
    concat_filename: str = "mite_concat.fasta"

    def run(self: Self) -> None:
        """Class entry point to run methods"""
        logger.info("BlastManager: Started")

        self.target_download.mkdir(exist_ok=True)
        self.target_blast.mkdir(exist_ok=True)

        if self.target_download.joinpath("MiteBlastDB.zip").exists():
            logger.warning("BlastManager: file 'MiteBlastDB' already exists - SKIP")
            return

        self.concat_fasta_files()
        self.generate_blast_db()

        logger.info("BlastManager: Completed")

    def concat_fasta_files(self: Self) -> None:
        """Concatenates individual FASTA files into a single one for BLAST processing."""
        with open(self.target_download.joinpath(self.concat_filename), "w") as outfile:
            for filename in self.src.iterdir():
                if filename.suffix == ".fasta" and filename != self.concat_filename:
                    with open(filename) as infile:
                        shutil.copyfileobj(infile, outfile)
                        outfile.write("\n")

    def generate_blast_db(self: Self) -> None:
        """Starts subprocess to generate a BLAST DB from the (downloaded) protein FASTA files"""
        with open(Path(__file__).parent.parent.joinpath("data/version.json")) as infile:
            version_data = json.load(infile)
            version = version_data.get("version_mite_data")

        command = [
            "makeblastdb",
            "-in",
            f"{self.target_download.joinpath(self.concat_filename)}",
            "-dbtype",
            "prot",
            "-out",
            f"{self.target_blast.joinpath("mite_blastfiles")}",
            "-title",
            f"MITE v{version} BLAST DB",
        ]
        subprocess.run(command, check=True)

        shutil.make_archive(
            base_name=str(self.target_download.joinpath("MiteBlastDB")),
            format="zip",
            root_dir=self.target_blast,
            base_dir=".",
        )
