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
    status = db.Column(db.String)
    retirement_reasons = db.Column(db.Text, nullable=True)
    comment = db.Column(db.Text, nullable=True)

    changelogs = db.relationship("ChangeLog", back_populates="entry")
    enzyme = db.relationship("Enzyme", uselist=False, back_populates="entry")


class Person(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    orcid = db.Column(db.String, unique=True)

    contributions = db.relationship(
        "ChangeLog", secondary="changelog_contributors", back_populates="contributors"
    )
    reviews = db.relationship(
        "ChangeLog", secondary="changelog_reviewers", back_populates="reviewers"
    )


class ChangeLog(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    entry_id = db.Column(db.String, db.ForeignKey("entry.accession"))
    entry = db.relationship("Entry", back_populates="changelogs")

    version = db.Column(db.String)
    date = db.Column(db.Date)
    comment = db.Column(db.Text)

    contributors = db.relationship(
        "Person", secondary="changelog_contributors", back_populates="contributions"
    )
    reviewers = db.relationship(
        "Person", secondary="changelog_reviewers", back_populates="reviews"
    )


# TODO: why not primary_key=True?
changelog_contributors = db.Table(
    "changelog_contributors",
    db.Column("person_id", db.Integer, db.ForeignKey("person.id")),
    db.Column("changelog_id", db.Integer, db.ForeignKey("change_log.id")),
)
# TODO: why not primary_key=True?
changelog_reviewers = db.Table(
    "changelog_reviewers",
    db.Column("person_id", db.Integer, db.ForeignKey("person.id")),
    db.Column("changelog_id", db.Integer, db.ForeignKey("change_log.id")),
)


class Cofactor(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    cofactor_name = db.Column(db.String, unique=True)
    cofactor_type = db.Column(db.String)

    enzymes = db.relationship(
        "Enzyme", secondary="enzyme_cofactors", back_populates="cofactors"
    )


class Reference(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    doi = db.Column(db.String, unique=True)

    enzymes = db.relationship(
        "Enzyme", secondary="enzyme_references", back_populates="references"
    )


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

    cofactors = db.relationship(
        "Cofactor", secondary="enzyme_cofactors", back_populates="enzymes"
    )
    references = db.relationship(
        "Reference", secondary="enzyme_references", back_populates="enzymes"
    )


enzyme_cofactors = db.Table(
    "enzyme_cofactors",
    db.Column(
        "cofactor_id", db.Integer, db.ForeignKey("cofactor.id"), primary_key=True
    ),
    db.Column("enzyme_id", db.Integer, db.ForeignKey("enzyme.id"), primary_key=True),
)

enzyme_references = db.Table(
    "enzyme_references",
    db.Column(
        "reference_id", db.Integer, db.ForeignKey("reference.id"), primary_key=True
    ),
    db.Column("enzyme_id", db.Integer, db.ForeignKey("enzyme.id"), primary_key=True),
)
