"""Generate auxiliary and essential files for webpage

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
import pickle
import shutil
from pathlib import Path

import pandas as pd
import requests
from Bio import Entrez, SeqIO
from pydantic import BaseModel
from rdkit.Chem import PandasTools, rdChemReactions

Entrez.email = "your_email@example.com"  # must be set but does not have to be real

logger = logging.getLogger("prep_data")


class Locations(BaseModel):
    """Base class to define the paths to storage directories

    Attributes:
        src: the location in which the downloaded mite json files are stored
        target: the location in which the auxiliary files are stored
        download: the location in which download-files are stored
    """

    src: Path = Path(__file__).parent.parent.joinpath("data/data")
    target: Path = Path(__file__).parent.parent.joinpath("data/")
    download: Path = Path(__file__).parent.parent.joinpath("data/download/")


class SummaryManager(Locations):
    """Prepare the data summary file

    Attributes:
        summary: dict to assemble data for later dump
        csv: dict to construct table for dump
        active_files: list of active mite entries (non-retired)
    """

    summary: dict = {"entries": {}}
    csv: dict = {
        "accession": [],
        "status": [],
        "name": [],
        "tailoring": [],
        "cofactors_organic": [],
        "cofactors_inorganic": [],
        "description": [],
        "reaction_description": [],
        "organism": [],
        "domain": [],
        "kingdom": [],
        "phylum": [],
        "class": [],
        "order": [],
        "family": [],
    }
    active_files: list = []

    def run(self) -> None:
        """Runs methods"""
        logger.info("SummaryManager: Started")

        if self.target.joinpath("summary.json").exists():
            logger.warning(
                f"SummaryManager: File {self.target.joinpath("summary.json")} already exists - skip file preparation"
            )
            return

        self.assemble_summary()
        self.dump_files()

        logger.info("SummaryManager: Completed")

    def assemble_summary(self):
        """Assemble the summary data"""
        for entry in self.src.iterdir():
            with open(entry) as infile:
                mite_data = json.load(infile)

            if mite_data.get("status") == "active":
                self.active_files.append(f"{mite_data["accession"]}.json")

            origin = self.get_organism(data=mite_data)

            self.summary["entries"][mite_data.get("accession")] = {
                "accession": mite_data["accession"],
                "status": '<i class="bi bi-check-circle-fill"></i>'
                if mite_data.get("status") == "active"
                else '<i class="bi bi-circle"></i>',
                "status_plain": mite_data.get("status"),
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
                "cofactors_organic": "N/A",
                "cofactors_inorganic": "N/A",
                "description": mite_data.get("enzyme", {}).get(
                    "description", "No description"
                ),
                "reaction_description": mite_data["reactions"][0].get(
                    "description", "No description"
                ),
                "organism": origin["organism"],
                "domain": origin["domain"],
                "kingdom": origin["kingdom"],
                "phylum": origin["phylum"],
                "class": origin["class"],
                "order": origin["order"],
                "family": origin["family"],
            }

            if (
                val := mite_data.get("enzyme", {})
                .get("cofactors", {})
                .get("organic", [])
            ):
                self.summary["entries"][mite_data.get("accession")][
                    "cofactors_organic"
                ] = "|".join(sorted(set(val)))

            if (
                val := mite_data.get("enzyme", {})
                .get("cofactors", {})
                .get("inorganic", [])
            ):
                self.summary["entries"][mite_data.get("accession")][
                    "cofactors_inorganic"
                ] = "|".join(sorted(set(val)))

        keys = sorted(self.summary.get("entries").keys())
        self.summary = {"entries": {key: self.summary["entries"][key] for key in keys}}

        for key, val in self.summary["entries"].items():
            self.csv["accession"].append(key)
            self.csv["status"].append(val["status_plain"])
            self.csv["name"].append(val["name"])
            self.csv["tailoring"].append(val["tailoring"])
            self.csv["cofactors_organic"].append(val["cofactors_organic"])
            self.csv["cofactors_inorganic"].append(val["cofactors_inorganic"])
            self.csv["description"].append(val["description"])
            self.csv["reaction_description"].append(val["reaction_description"])
            self.csv["organism"].append(val["organism"])
            self.csv["domain"].append(val["domain"])
            self.csv["kingdom"].append(val["kingdom"])
            self.csv["phylum"].append(val["phylum"])
            self.csv["class"].append(val["class"])
            self.csv["order"].append(val["order"])
            self.csv["family"].append(val["family"])

    @staticmethod
    def get_organism(data: dict) -> dict:
        """Download the organism identifier

        Arguments:
            data: mite entry

        Returns:
            A dict of organism info
        """
        origin = {
            "organism": "Not found",
            "domain": "Not found",
            "kingdom": "Not found",
            "phylum": "Not found",
            "class": "Not found",
            "order": "Not found",
            "family": "Not found",
        }

        if acc := data["enzyme"]["databaseIds"].get("uniprot"):
            if (
                response := requests.get(
                    f"https://rest.uniprot.org/uniprotkb/{acc}.json"
                )
            ).status_code == 200 or (
                response := requests.get(f"https://rest.uniprot.org/uniparc/{acc}.json")
            ).status_code == 200:
                content = response.json()
                record = content.get("organism", {}).get("lineage", None)
                try:
                    origin["domain"] = record[0] or "Not found"
                    origin["kingdom"] = record[1] or "Not found"
                    origin["phylum"] = record[2] or "Not found"
                    origin["class"] = record[3] or "Not found"
                    origin["order"] = record[4] or "Not found"
                    origin["family"] = record[5] or "Not found"
                    origin["organism"] = (
                        content.get("organism", {}).get("scientificName") or "Not found"
                    )
                except Exception as e:
                    logger.warning(
                        f"Did not find organism information for '{acc}' ({data["accession"]}): {e!s}"
                    )
                return origin

        if acc := data["enzyme"]["databaseIds"].get("genpept"):
            try:
                handle = Entrez.efetch(
                    db="protein", id=acc, rettype="gb", retmode="text"
                )
                record = SeqIO.read(handle, "genbank")
                origin["domain"] = record.annotations.get("taxonomy")[0] or "Not found"
                origin["kingdom"] = record.annotations.get("taxonomy")[1] or "Not found"
                origin["phylum"] = record.annotations.get("taxonomy")[2] or "Not found"
                origin["class"] = record.annotations.get("taxonomy")[3] or "Not found"
                origin["order"] = record.annotations.get("taxonomy")[4] or "Not found"
                origin["family"] = record.annotations.get("taxonomy")[5] or "Not found"
                origin["organism"] = record.annotations.get("organism") or "Not found"
                handle.close()
            except Exception as e:
                logger.warning(
                    f"Did not find organism information for '{acc}' ({data["accession"]}): {e!s}"
                )
            return origin

        return origin

    def dump_files(self):
        """Dump files to disk"""
        with open(self.target.joinpath("summary.json"), "w") as outfile:
            outfile.write(json.dumps(self.summary, indent=2, ensure_ascii=False))

        df = pd.DataFrame(self.csv)
        df.to_csv(self.target.joinpath("summary.csv"), index=False)

        self.download.mkdir(parents=True, exist_ok=True)

        self.download.joinpath("temp_dir").mkdir(parents=True)
        for filename in self.active_files:
            shutil.copy(
                src=self.src.joinpath(filename),
                dst=self.download.joinpath(f"temp_dir/{filename}"),
            )
        shutil.make_archive(
            base_name=str(self.download.joinpath("MITE_all_active_entries").resolve()),
            format="zip",
            root_dir=self.download.joinpath("temp_dir"),
            base_dir=".",
        )
        shutil.rmtree(self.download.joinpath("temp_dir"))


class AuxFileManager(Locations):
    """Prepare auxiliary files for website

    Only entries with "active" flag are used to compile auxiliary files

    Attributes:
        smiles: dict to SMILES be exported as csv file
        smarts: dict of reaction SMARTS to be exported as csv file
        pickle_substrates: list with pre-calculated fingerprints for substructure search
        pickle_products: list with pre-calculated fingerprints for substructure search
        pickle_smartsfps: dict with pre-calucated reaction smarts fingerprints for search
    """

    smiles: dict = {"mite_id": [], "substrates": [], "products": []}
    smarts: dict = {
        "mite_id": [],
        "reactionsmarts": [],
    }
    pickle_substrates: list = []
    pickle_products: list = []
    pickle_smartsfps: dict = {
        "mite_id": [],
        "reactionsmarts": [],
        "reaction_fps": [],
    }

    def run(self) -> None:
        """Call methods for preparation of auxiliary files"""
        logger.info("AuxFileManager: Started")

        if self.download.joinpath("dump_smiles.csv").exists():
            logger.warning(
                f"AuxFileManager: File {self.download.joinpath("dump_smiles.csv")} already exists - skip file preparation"
            )
            return

        self.prepare_files()
        self.dump_files()

        logger.info("AuxFileManager: Completed")

    def prepare_files(self) -> None:
        """Prepare the auxiliary files derived from mite entries

        Only "active" entries are dumped; others are skipped
        """
        self.download.mkdir(parents=True, exist_ok=True)

        for entry in self.src.iterdir():
            with open(entry) as infile:
                mite_data = json.load(infile)

            if mite_data["status"] != "active":
                continue

            self.prepare_smiles(mite_data)
            self.prepare_smarts(mite_data)

        self.prepare_pickled_smiles()
        self.prepare_pickled_smarts()

    def prepare_smiles(self, data: dict) -> None:
        """Create a table of SMILES strings contained in MITE entries

        Arguments:
            data: a dict derived from a mite json file
        """
        for readctionid, reaction in enumerate(data["reactions"], 1):
            for exampleid, example in enumerate(reaction["reactions"], 1):
                self.smiles["mite_id"].append(
                    f"{data['accession']}.reaction{readctionid}.example{exampleid}"
                )
                self.smiles["substrates"].append(f"{example['substrate']}")
                self.smiles["products"].append(f"{'.'.join(example['products'])}")

    def prepare_smarts(self, data: dict) -> None:
        """Create a table of reaction SMARTS strings contained in MITE entries

        Arguments:
            data: a dict derived from a mite json file
        """
        for readctionid, reaction in enumerate(data["reactions"], 1):
            self.smarts["mite_id"].append(f'{data['accession']}.reaction{readctionid}')
            self.smarts["reactionsmarts"].append(f"{reaction["reactionSMARTS"]}")

    def prepare_pickled_smiles(self) -> None:
        """Create a pickle file with pre-calculated SMILES fingerprints"""
        df = pd.DataFrame(self.smiles)

        PandasTools.AddMoleculeColumnToFrame(
            df,
            smilesCol="substrates",
            molCol="ROMol_substrates",
            includeFingerprints=True,
        )
        PandasTools.AddMoleculeColumnToFrame(
            df, smilesCol="products", molCol="ROMol_products", includeFingerprints=True
        )

        self.pickle_substrates = list(df["ROMol_substrates"])
        self.pickle_products = list(df["ROMol_products"])

    def prepare_pickled_smarts(self) -> None:
        """Create a pickle file of a dict with pre-calculated reaction SMARTS fingerprints"""
        self.pickle_smartsfps["mite_id"] = self.smarts.get("mite_id")
        self.pickle_smartsfps["reactionsmarts"] = self.smarts.get("reactionsmarts")

        for smarts in self.smarts.get("reactionsmarts"):
            self.pickle_smartsfps["reaction_fps"].append(
                rdChemReactions.CreateStructuralFingerprintForReaction(
                    rdChemReactions.ReactionFromSmarts(smarts)
                )
            )

    def dump_files(self) -> None:
        """Dump the assembled files"""
        df_smiles = pd.DataFrame(self.smiles)
        df_smiles.to_csv(path_or_buf=self.download.joinpath("dump_smiles.csv"))

        df_smarts = pd.DataFrame(self.smarts)
        df_smarts.to_csv(path_or_buf=self.download.joinpath("dump_smarts.csv"))

        with open(self.target.joinpath("substrate_list.pickle"), "wb") as outfile:
            pickle.dump(obj=self.pickle_substrates, file=outfile)

        with open(self.target.joinpath("product_list.pickle"), "wb") as outfile:
            pickle.dump(obj=self.pickle_products, file=outfile)

        with open(self.target.joinpath("reaction_fps.pickle"), "wb") as outfile:
            pickle.dump(obj=self.pickle_smartsfps, file=outfile)
