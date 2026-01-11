import pickle
from copy import deepcopy
from pathlib import Path
from typing import Literal

import pandas as pd
from pydantic import BaseModel, Field
from rdkit.Chem import (
    DataStructs,
    PandasTools,
    rdChemReactions,
)
from rdkit.DataStructs import FingerprintSimilarity

from app.core.config import settings


def load_pickle(p: Path):
    with open(p, "rb") as fin:
        return pickle.load(fin)


class ReactionProvider:
    """Lazy-load reaction fingerprints"""

    _fp_reaction = None

    @classmethod
    def get(cls):
        if cls._fp_reaction is None:
            cls._fp_reaction = load_pickle(
                settings.data_dir.joinpath("reaction_fps.pickle")
            )
        return cls._fp_reaction


class ReactionSearch(BaseModel):
    """Organizes querying for a reaction (SMARTS)"""

    reaction: str = Field(..., description="Reaction query")
    sim_score: float = Field(..., description="Similarity score threshold")
    mode: Literal["rd_difference", "rd_structural"]

    def search_reaction(self) -> dict:
        """Query dataset for a specific reaction"""

        reaction_dict = deepcopy(ReactionProvider.get())

        if self.mode == "rd_structural":
            reaction_fp = rdChemReactions.CreateStructuralFingerprintForReaction(
                rdChemReactions.ReactionFromSmarts(self.reaction)
            )
            similarities = [
                FingerprintSimilarity(reaction_fp, fp)
                for fp in reaction_dict["reaction_fps"]
            ]
            reaction_dict["similarities"] = similarities
        else:
            reaction_fp = rdChemReactions.CreateDifferenceFingerprintForReaction(
                rdChemReactions.ReactionFromSmarts(self.reaction)
            )
            similarities = [
                DataStructs.TanimotoSimilarity(reaction_fp, fp)
                for fp in reaction_dict["diff_reaction_pfs"]
            ]
            reaction_dict["similarities"] = similarities

        summary = {}
        df = pd.DataFrame(reaction_dict)
        for _, row in df.iterrows():
            key = row["mite_id"].split(".")[0]
            if row["similarities"] >= self.sim_score:
                if key in summary:
                    summary[key]["reaction"].append(
                        row["mite_id"].split(".")[1].replace("reaction", "rx")
                    )
                    summary[key]["sim_score"].append(round(row["similarities"], 2))
                else:
                    summary[key] = {
                        "reaction": [
                            row["mite_id"].split(".")[1].replace("reaction", "rx")
                        ],
                        "sim_score": [round(row["similarities"], 2)],
                    }

        return summary
