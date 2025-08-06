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


changelog_contributors = db.Table(
    "changelog_contributors",
    db.Column("person_id", db.Integer, db.ForeignKey("person.id")),
    db.Column("changelog_id", db.Integer, db.ForeignKey("change_log.id")),
)

changelog_reviewers = db.Table(
    "changelog_reviewers",
    db.Column("person_id", db.Integer, db.ForeignKey("person.id")),
    db.Column("changelog_id", db.Integer, db.ForeignKey("change_log.id")),
)
