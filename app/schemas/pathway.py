import base64
import csv
import io
import re
from typing import Annotated, Any, ClassVar

from pydantic import (
    BaseModel,
    Field,
    PositiveInt,
    computed_field,
    field_validator,
    model_validator,
)
from rdkit.Chem import Draw, MolFromSmiles, MolToSmiles, rdChemReactions

from app.core.config import settings
from app.services.file_handling import load_json

MITE_RE = re.compile(r"^MITE(\d{7})$")


class PathwayParams(BaseModel):
    """Holds in silico pathway processing parameters

    Attributes:
        substrate: An enzyme substrate (structure) in SMILES format
    """

    substrate: str
    mite_id: list[Annotated[str, Field(pattern=MITE_RE)]]
    reaction_nr: list[PositiveInt]

    @field_validator("substrate")
    @classmethod
    def validate_substrate(cls, v: str):
        if not MolFromSmiles(v):
            raise ValueError("Substrate is not a valid SMILES string.")
        return v

    @model_validator(mode="after")
    def validate_alignment(self):
        if len(self.mite_id) == 0:
            raise ValueError("At least one MITE ID must be provided")

        if len(self.reaction_nr) == 0:
            raise ValueError("At least one reaction number must be provided")

        if len(self.mite_id) != len(self.reaction_nr):
            raise ValueError("mite_id and reaction_nr must have the same length.")
        return self

    @computed_field
    @property
    def reactions(self) -> list[tuple[str, int]]:
        """Create mapping MITE ID -> reaction #"""
        return list(zip(self.mite_id, self.reaction_nr, strict=False))

    def dump_params(self) -> dict:
        """Format params for params dumping and Jinja2 processing"""
        params_dict = self.model_dump()

        for k, v in params_dict.items():
            if not isinstance(v, list):
                params_dict[k] = [v]

        return params_dict


class RunPathway(BaseModel):
    svgs: list = []
    report: dict = {"substrate": [], "MITE_acc": [], "enzyme_name": [], "product": []}

    @staticmethod
    def generate_svg(mol: Any) -> str:
        smiles = MolFromSmiles(mol)
        for atom in smiles.GetAtoms():
            atom.SetAtomMapNum(0)
        smiles = Draw.rdMolDraw2D.PrepareMolForDrawing(smiles)
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(-1, -1)
        dopts = drawer.drawOptions()
        dopts.clearBackground = False
        drawer.DrawMolecule(smiles)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return base64.b64encode(svg.encode("utf-8")).decode("utf-8")

    def get_pathway(self, p: PathwayParams) -> None:
        """Create pathway by applying reactions to substrate step-by-step"""

        substrate = p.substrate
        self.svgs.append(self.generate_svg(mol=substrate))

        for step in p.reactions:
            data = load_json(settings.data_dir.joinpath(f"data/{step[0]}.json"))

            index = int(step[1]) - 1
            if index > len(data["reactions"]):
                raise IndexError(f"Reaction does not exist in {step[0]}: {step[1]}")

            rd_substrate = MolFromSmiles(substrate)
            rd_reaction = rdChemReactions.ReactionFromSmarts(
                data["reactions"][index]["reactionSMARTS"]
            )
            products = rd_reaction.RunReactants([rd_substrate])
            products = {MolToSmiles(product[0]) for product in products}
            products = list(products)
            products.sort()

            self.report["substrate"].append(substrate)
            self.report["MITE_acc"].append(f"{step[0]}:{step[1]}")
            self.report["enzyme_name"].append(data["enzyme"]["name"])

            try:
                if len(products) != 0:
                    self.report["product"].append(products[0])
                    self.svgs.append(self.generate_svg(products[0]))
                else:
                    self.report["product"].append("")
                    break
            except Exception as e:
                raise ValueError(f"Error during reaction {step[0]}: {step[1]}:") from e

            substrate = products[0]

    def stream_csv(self) -> io.StringIO:
        if not self.report:
            return

        print(self.report)

        buffer = io.StringIO()
        writer = csv.writer(buffer)

        writer.writerow(self.report.keys())
        yield buffer.getvalue()
        buffer.seek(0)
        buffer.truncate(0)

        for row in zip(*self.report.values(), strict=False):
            writer.writerow(row)
            yield buffer.getvalue()
            buffer.seek(0)
            buffer.truncate(0)
