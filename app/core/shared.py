import json

from mite_schema import SchemaManager

from app.core.config import settings
from app.services.file_handling import load_json


def load_table_head():
    """Load headers for overview table"""
    return [
        ("accession", "MITE Accession"),
        ("status", "Status"),
        ("name", "Enzyme Name"),
        ("tailoring", "Reaction Tailoring Term"),
        ("description", "Enzyme Description"),
        ("organism", "Taxonomy Organism"),
        ("family", "Taxonomy Family"),
        ("reaction_description", "Reaction Description"),
    ]


def load_active():
    """Load active (=non-retired) entries"""
    return load_json(settings.data_dir.joinpath("active.json"))


def load_retired():
    """Load inactive (=retired) entries"""
    return load_json(settings.data_dir.joinpath("retired.json"))


def load_form_vals():
    """Load schema definitions"""
    with open(SchemaManager().entry) as f:
        schema = json.load(f)
        return {
            "evidence": schema["$defs"]["evidence"]["enum"],
            "tailoring": schema["$defs"]["tailoringFunction"]["enum"],
            "inorganic": schema["$defs"]["inorganic"]["enum"],
            "organic": schema["$defs"]["organic"]["enum"],
        }
