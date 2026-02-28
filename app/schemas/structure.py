import pickle
from pathlib import Path
from typing import Literal

import pandas as pd
from pydantic import BaseModel, Field
from rdkit.Chem import MolFromSmarts, MolFromSmiles, PandasTools

from app.core.config import settings


def load_pickle(p: Path):
    with open(p, "rb") as fin:
        return pickle.load(fin)


class SubstratesProvider:
    """Lazy-load substrate fingerprints"""

    _fp_substrates = None

    @classmethod
    def get(cls):
        if cls._fp_substrates is None:
            cls._fp_substrates = load_pickle(
                settings.data_dir.joinpath("substrate_list.pickle")
            )
        return cls._fp_substrates


class ProductsProvider:
    """Lazy-load products fingerprints"""

    _fp_products = None

    @classmethod
    def get(cls):
        if cls._fp_products is None:
            cls._fp_products = load_pickle(
                settings.data_dir.joinpath("product_list.pickle")
            )
        return cls._fp_products


class SmilesProvider:
    """Lazy-load SMILES"""

    _df = None

    @classmethod
    def get(cls):
        """Merges columns back into dataframe"""
        if cls._df is None:
            df = pd.read_csv(settings.data_dir.joinpath("download/dump_smiles.csv"))
            df["ROMol_substrates"] = SubstratesProvider.get()
            df["ROMol_products"] = ProductsProvider.get()
            cls._df = df
        return cls._df


class StructureSearch(BaseModel):
    """Organizes querying for a smiles/smarts (sub)structure"""

    structure: str = Field(..., description="Substructure query")
    s_type: Literal["smiles", "smarts"]

    def search_structure(self) -> dict:
        """Query dataset for a specific (sub)structure"""
        df = SmilesProvider.get()

        if self.s_type == "smiles":
            mol_query = MolFromSmiles(self.structure)
        else:
            mol_query = MolFromSmarts(self.structure)

        df_substrate_match = df[df["ROMol_substrates"] >= mol_query]
        df_product_match = df[df["ROMol_products"] >= mol_query]

        matches = set()
        matches.update(set(df_substrate_match["mite_id"]))
        matches.update(set(df_product_match["mite_id"]))

        matches = sorted(matches)

        summary = {}
        for match in matches:
            reaction = match.split(".")[1].replace("reaction", "rx")
            example = match.split(".")[2].replace("example", "ex")
            key = match.split(".")[0]

            if key in summary:
                summary[key]["reaction"].append(reaction)
                summary[key]["example"].append(f"{reaction}-{example}")
            else:
                summary[key] = {
                    "reaction": [reaction],
                    "example": [f"{reaction}-{example}"],
                }

        return summary
