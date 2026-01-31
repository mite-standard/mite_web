from typing import Any

from sqlalchemy import ColumnElement, and_, inspect, or_
from sqlalchemy.orm import Session

from app.db.models import Entry, Enzyme

operators = {
    "equal": lambda col, val: col == val,
    "not_equal": lambda col, val: col != val,
    "contains": lambda col, val: col.ilike(f"%{val}%"),
    "not_contains": lambda col, val: ~col.ilike(f"%{val}%"),
    "is_null": lambda col, _: col.is_(None),
    "is_not_null": lambda col, _: col.is_not(None),
}

field_map = {
    "orcids",
    "references",
    "evidences",
    "tailoring",
    "enzyme.name",
    "enzyme.enzyme_description",
    "enzyme.uniprot_id",
    "enzyme.genpept_id",
    "enzyme.mibig_id",
    "enzyme.wikidata_id",
    "enzyme.has_auxenzymes",
    "enzyme.organism_id",
    "enzyme.domain_id",
    "enzyme.kingdom_id",
    "enzyme.phylum_id",
    "enzyme.class_id",
    "enzyme.order_id",
    "enzyme.family_id",
    "enzyme.cofactors",
    "reactions.description",
    "reactions.reaction_smarts",
    "reactions.rhea_id",
    "reactions.ec_id",
    "reactions.example_reactions.smiles_substrate",
    "reactions.example_reactions.is_intermediate",
    "reactions.example_reactions.products.smiles_product",
}


def enzyme_exists(accession: str, db: Session) -> str | None:
    """Query genpept/uniprot fields for accession, return MITE accession"""
    match = (
        db.query(Entry)
        .join(Entry.enzyme)
        .filter(or_(Enzyme.uniprot_id == accession, Enzyme.genpept_id == accession))
        .first()
    )
    if match:
        return match.accession


def search_db(rules: dict, db: Session) -> set[str]:
    """Queries the DB

    Args:
        rules: the raw form from QueryBuilder JSON rules
        db: the database session

    Returns:
        A set of MITE accession IDs for filtering
    """
    filters = parse_rules_to_filters(rules, Entry)
    query = db.query(Entry).filter(filters)
    return {e.accession for e in query}


def parse_rules_to_filters(rules: dict, base_model: Any) -> ColumnElement:
    """Convert QueryBuilder JSON rules into a SQLAlchemy-compatible expression

    Note: all expressions are chained together with either 'AND' or 'OR'
    """

    expressions = []
    for rule in rules.get("rules", []):
        field_path = rule["field"] if rule["field"] in field_map else None
        if not field_path:
            continue

        filter_expr = build_filter_from_path(
            model=base_model,
            path_parts=field_path.split("."),
            operator=rule["operator"],
            value=rule.get("value"),
        )
        expressions.append(filter_expr)

    if rules.get("condition", "AND").upper() == "OR":
        return or_(*expressions)
    else:
        return and_(*expressions)


def build_filter_from_path(model: Any, path_parts: list, operator: str, value: str):
    """Build a SQLAlchemy filter from a dotted path."""
    mapper = inspect(model)

    if path_parts[0] in mapper.relationships:
        rel = mapper.relationships[path_parts[0]]
        related_model = rel.entity.class_

        nested_filter = build_filter_from_path(
            related_model, path_parts[1:], operator, value
        )

        if rel.uselist:
            return rel.class_attribute.any(nested_filter)
        else:
            return rel.class_attribute.has(nested_filter)
    elif path_parts[0] in mapper.columns:
        col = mapper.columns[path_parts[0]]
        return operators[operator](col, value)
    else:
        raise ValueError(f"Invalid path segment: {path_parts[0]}")
