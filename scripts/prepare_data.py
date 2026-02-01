import json
import logging
import os
import shutil
import sys
import time
from importlib import metadata
from pathlib import Path
from typing import Any

import requests
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
if not logger.handlers:
    logger.addHandler(handler)


def get_version() -> str:
    try:
        return metadata.version("mite_web")
    except metadata.PackageNotFoundError:
        return "0+unknown"


class RecordManager(BaseModel):
    """Download data and prepare for use by mite_web

    Attributes:
        r_data: Record for mite_data zenodo
        r_extras: Record for mite_web_extras zenodo
        loc: data location
        session: a Request Session to keep code DRY
    """

    r_data: str
    r_extras: str
    loc: Path = Field(
        default_factory=lambda: Path(__file__).parent.parent.joinpath("data")
    )
    session: Any = Field(default_factory=requests.Session)

    def __init__(self, **data):
        super().__init__(**data)
        self.session.headers.setdefault(
            "User-Agent",
            f"mite_web/{get_version()} (https://github.com/mite-standard/mite_web)",
        )

    def run(self) -> None:
        """Call methods for downloading and moving data"""

        logger.info("RecordManager: Started")
        self.loc.mkdir()

        self.download_mite_data()
        self.organize_mite_data()

        self.download_mite_web_extras()
        self.organize_mite_web_extras()

        self.generate_summary()
        self.diff_active_retired()

        logger.info("RecordManager: Completed")

    def download_mite_data(self) -> None:
        """Download and arrange mite_data record from Zenodo

        Raises:
            RuntimeError: Could not download files
        """
        logger.info("RecordManager: Start mite_data download")

        rsps_meta = self.session.get(
            f"https://zenodo.org/api/records/{self.r_data}", timeout=12.1
        )
        if rsps_meta.status_code != 200:
            raise RuntimeError(
                f"Error fetching 'mite_data' metadata: {rsps_meta.status_code}. Reason: {rsps_meta.reason}"
            )
        record_metadata = rsps_meta.json()
        version = record_metadata["metadata"]["version"]
        files_url = record_metadata["files"][0]["links"]["self"]

        time.sleep(1)

        rsps_data = self.session.get(files_url, timeout=60.1)
        if rsps_data.status_code != 200:
            raise RuntimeError(
                f"Error downloading 'mite_data' record: {rsps_data.status_code}. Reason: {rsps_data.reason}"
            )

        with open(self.loc.joinpath("mite_data.zip"), "wb") as f:
            f.write(rsps_data.content)

        with open(self.loc.joinpath("version_mite_data.json"), "w") as f:
            f.write(json.dumps({"version_mite_data": f"{version}"}))

        logger.info("RecordManager: Completed mite_data download")

    def download_mite_web_extras(self) -> None:
        """Download and arrange mite_web_extras record from Zenodo

        Raises:
            RuntimeError: Could not download files
        """
        logger.info("RecordManager: Start mite_web_extras download")

        rsps_meta = self.session.get(
            f"https://zenodo.org/api/records/{self.r_extras}", timeout=12.1
        )
        if rsps_meta.status_code != 200:
            raise RuntimeError(
                f"Error fetching 'mite_web_extras' metadata: {rsps_meta.status_code}. Reason: {rsps_meta.reason}"
            )
        record_metadata = rsps_meta.json()
        files_url = record_metadata["files"][0]["links"]["self"]

        time.sleep(1)

        rsps_data = self.session.get(files_url, timeout=60.1)
        if rsps_data.status_code != 200:
            raise RuntimeError(
                f"Error downloading 'mite_web_extras' record: {rsps_data.status_code}. Reason: {rsps_data.reason}"
            )

        with open(self.loc.joinpath("mite_web_extras.zip"), "wb") as f:
            f.write(rsps_data.content)

        logger.info("RecordManager: Completed mite_web_extras download")

    def organize_mite_data(self) -> None:
        """Unpacks mite_data, moves to correct locations, cleans up

        Raises:
            NotADirectoryError: directory not unzipped in expected location
            RuntimeError: Could not determine data location in downloaded folder
        """
        logger.info("RecordManager: Organizing mite_data download")

        src = self.loc.joinpath("mite_data.zip")
        trgt = self.loc.joinpath("mite_data_unzip")
        download = self.loc.joinpath("download")

        download.mkdir(exist_ok=True)

        shutil.unpack_archive(filename=src, extract_dir=trgt, format="zip")
        if not trgt.exists():
            raise NotADirectoryError(f"Could not find the unzipped directory {trgt}.")

        matching_dirs = list(trgt.glob("mite-standard-mite_data-*"))
        if not matching_dirs:
            raise RuntimeError(
                f"Could not determine data storage location in downloaded directory."
            )
        sdir = matching_dirs[0]

        shutil.move(
            src=trgt.joinpath(sdir).joinpath("mite_data/data").resolve(),
            dst=self.loc.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir).joinpath("mite_data/fasta").resolve(),
            dst=self.loc.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir)
            .joinpath("mite_data/metadata/dump_smarts.csv")
            .resolve(),
            dst=download.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir)
            .joinpath("mite_data/metadata/dump_smiles.csv")
            .resolve(),
            dst=download.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir)
            .joinpath("mite_data/metadata/metadata_general.json")
            .resolve(),
            dst=self.loc.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir)
            .joinpath("mite_data/metadata/metadata_mibig.json")
            .resolve(),
            dst=self.loc.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir)
            .joinpath("mite_data/metadata/product_list.pickle")
            .resolve(),
            dst=self.loc.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir)
            .joinpath("mite_data/metadata/reaction_fps.pickle")
            .resolve(),
            dst=self.loc.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir)
            .joinpath("mite_data/metadata/substrate_list.pickle")
            .resolve(),
            dst=self.loc.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir)
            .joinpath("mite_data/metadata/summary.csv")
            .resolve(),
            dst=self.loc.resolve(),
        )

        os.remove(src)
        shutil.rmtree(trgt)

        logger.info("RecordManager: Completed organizing mite_data download")

    def organize_mite_web_extras(self) -> None:
        """Unpacks mite_web_extras, moves to correct locations, cleans up

        Raises:
            NotADirectoryError: directory not unzipped in expected location
            RuntimeError: Could not determine data location in downloaded folder
        """
        logger.info("RecordManager: Organizing mite_web_extras download")

        src = self.loc.joinpath("mite_web_extras.zip")
        trgt = self.loc.joinpath("mite_web_extras_unzip")
        download = self.loc.joinpath("download")

        download.mkdir(exist_ok=True)

        shutil.unpack_archive(filename=src, extract_dir=trgt, format="zip")
        if not trgt.exists():
            raise NotADirectoryError(f"Could not find the unzipped directory {trgt}.")

        matching_dirs = list(trgt.glob("mite-standard-mite_web_extras-*"))
        if not matching_dirs:
            raise RuntimeError(
                f"Could not determine data storage location in downloaded directory."
            )
        sdir = matching_dirs[0]

        shutil.move(
            src=trgt.joinpath(sdir).joinpath("data/blastlib").resolve(),
            dst=self.loc.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir).joinpath("data/html").resolve(),
            dst=self.loc.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir).joinpath("data/img").resolve(),
            dst=self.loc.resolve(),
        )
        shutil.move(
            src=trgt.joinpath(sdir).joinpath("data/mite_concat.fasta").resolve(),
            dst=download.resolve(),
        )

        with open(self.loc.joinpath("version_mite_data.json")) as infile:
            v_mite_data = json.load(infile)

        with open(trgt.joinpath(sdir).joinpath("data/version.json")) as infile:
            v_mite_web_extr = json.load(infile)

        if v_mite_data["version_mite_data"] != v_mite_web_extr["version_mite_data"]:
            raise RuntimeError(
                f"Version mite_data {v_mite_data["version_mite_data"]} does not match version mite_web_extras {v_mite_web_extr["version_mite_data"]}. \n"
                f"Was mite_web_extras updated together with mite_data?"
            )

        os.remove(src)
        shutil.rmtree(trgt)

        logger.info("RecordManager: Completed organization of mite_web_extras download")

    def generate_summary(self) -> None:
        """Creates summary file from metadata_general.json"""
        logger.info("RecordManager: Starting creation of summary")
        summary = {"entries": {}}

        with open(self.loc.joinpath("metadata_general.json")) as infile:
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

        with open(self.loc.joinpath("summary.json"), "w") as f:
            f.write(json.dumps(summary))
        logger.info("RecordManager: Completed creation of summary")

    def diff_active_retired(self) -> None:
        """Populate separate lists for active and retired entries"""
        logger.info("RecordManager: Starting creation of active/retired")
        with open(self.loc.joinpath("summary.json")) as infile:
            summary = json.load(infile)
            active = {
                key: val
                for key, val in summary["entries"].items()
                if val.get("status_plain") == "active"
            }
            retired = {
                key: val
                for key, val in summary["entries"].items()
                if val.get("status_plain") == "retired"
            }
            with open(self.loc.joinpath("active.json"), "w") as f:
                f.write(json.dumps(active))
            with open(self.loc.joinpath("retired.json"), "w") as f:
                f.write(json.dumps(retired))
        logger.info("RecordManager: Completed creation of active/retired")


def main(data: str, extras: str):
    """Prepares mite_data for use in mite_web

    Args:
        data: the Zenodo mite_data record
        extras: the Zenodo mite_web_extras record

    Raises:
        RuntimeError: no Zenodo records specified
    """
    logger.info("Starting data prep.")

    if not data or not extras:
        raise RuntimeError(
            "Run as 'python prepare_data.py <mite_data record> <mite_web_extras record>'."
        )

    download_manager = RecordManager(r_data=data, r_extras=extras)
    download_manager.run()

    logger.info("Completed data prep.")


if __name__ == "__main__":
    main(data=sys.argv[1], extras=sys.argv[2])
