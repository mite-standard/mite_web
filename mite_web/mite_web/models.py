"""Creates schema for PostgreSQL database

For schema, see also:
https://github.com/mite-standard/mite_schema/blob/main/mite_schema/schema/entry.json

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

from mite_web.config.extensions import db


class Entry(db.Model):
    accession = db.Column(db.String, primary_key=True)

    orcids = db.Column(db.Text)
    references = db.Column(db.Text)

    enzyme = db.relationship("Enzyme", uselist=False, back_populates="entry")
    reactions = db.relationship("Reaction", back_populates="entry")

class Enzyme(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    entry_id = db.Column(db.String, db.ForeignKey("entry.accession"))
    entry = db.relationship("Entry", back_populates="enzyme")

    name = db.Column(db.Text)
    enzyme_description = db.Column(db.Text, nullable=True)
    uniprot_id = db.Column(db.String, nullable=True, index=True)
    genpept_id = db.Column(db.String, nullable=True, index=True)
    mibig_id = db.Column(db.String, nullable=True, index=True)
    wikidata_id = db.Column(db.String, nullable=True, index=True)
    has_auxenzymes = db.Column(db.Boolean, default=False)
    organism_id = db.Column(db.String, nullable=True, index=True)
    domain_id = db.Column(db.String, nullable=True, index=True)
    kingdom_id = db.Column(db.String, nullable=True, index=True)
    phylum_id = db.Column(db.String, nullable=True, index=True)
    class_id = db.Column(db.String, nullable=True, index=True)
    order_id = db.Column(db.String, nullable=True, index=True)
    family_id = db.Column(db.String, nullable=True, index=True)
    cofactors = db.Column(db.String, nullable=True, index=True)


class Tailoring(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    tailoring = db.Column(db.String, unique=True)

    reaction = db.relationship(
        "Reaction", secondary="reaction_tailoring", back_populates="tailorings"
    )


class Evidence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    evidence = db.Column(db.Text, unique=True)

    reaction = db.relationship(
        "Reaction", secondary="reaction_evidence", back_populates="evidences"
    )


class Product(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    example_reaction_id = db.Column(db.Integer, db.ForeignKey("example_reaction.id"))

    smiles_product = db.Column(db.Text)

    example_reaction = db.relationship("ExampleReaction", back_populates="products")


class ExampleReaction(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    reaction_id = db.Column(db.Integer, db.ForeignKey("reaction.id"))
    reaction = db.relationship("Reaction", back_populates="example_reactions")

    smiles_substrate = db.Column(db.Text)
    is_intermediate = db.Column(db.Boolean)

    products = db.relationship(
        "Product", back_populates="example_reaction", cascade="all, delete-orphan"
    )


class Reaction(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    entry_id = db.Column(db.String, db.ForeignKey("entry.accession"))
    entry = db.relationship("Entry", back_populates="reactions")

    description = db.Column(db.Text, nullable=True)
    reaction_smarts = db.Column(db.Text)
    rhea_id = db.Column(db.Integer, nullable=True, index=True)
    ec_id = db.Column(db.String, nullable=True, index=True)

    tailorings = db.relationship(
        "Tailoring", secondary="reaction_tailoring", back_populates="reaction"
    )
    evidences = db.relationship(
        "Evidence", secondary="reaction_evidence", back_populates="reaction"
    )
    example_reactions = db.relationship(
        "ExampleReaction", back_populates="reaction", cascade="all, delete-orphan"
    )


reaction_tailoring = db.Table(
    "reaction_tailoring",
    db.Column(
        "tailoring_id", db.Integer, db.ForeignKey("tailoring.id"), primary_key=True
    ),
    db.Column(
        "reaction_id", db.Integer, db.ForeignKey("reaction.id"), primary_key=True
    ),
)

reaction_evidence = db.Table(
    "reaction_evidence",
    db.Column(
        "evidence_id", db.Integer, db.ForeignKey("evidence.id"), primary_key=True
    ),
    db.Column(
        "reaction_id", db.Integer, db.ForeignKey("reaction.id"), primary_key=True
    ),
)

