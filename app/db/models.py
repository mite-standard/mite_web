from sqlalchemy import Boolean, Column, ForeignKey, Integer, String, Text
from sqlalchemy.orm import relationship

from app.db.database import Base


class Entry(Base):
    """Parent table describing individual entries"""

    __tablename__ = "entry"

    accession = Column(String, primary_key=True)
    orcids = Column(Text)
    references = Column(Text)
    evidences = Column(Text)
    tailoring = Column(Text)

    enzyme = relationship("Enzyme", uselist=False, back_populates="entry")
    reactions = relationship("Reaction", back_populates="entry")


class Enzyme(Base):
    """Child of Entry (one enzyme per entry)"""

    __tablename__ = "enzyme"

    id = Column(Integer, primary_key=True)
    entry_id = Column(String, ForeignKey("entry.accession"))
    entry = relationship("Entry", back_populates="enzyme")

    name = Column(Text)
    enzyme_description = Column(Text, nullable=True)
    uniprot_id = Column(String, nullable=True, index=True)
    genpept_id = Column(String, nullable=True, index=True)
    mibig_id = Column(String, nullable=True, index=True)
    wikidata_id = Column(String, nullable=True, index=True)
    has_auxenzymes = Column(Boolean, default=False)
    organism_id = Column(String, nullable=True, index=True)
    domain_id = Column(String, nullable=True, index=True)
    kingdom_id = Column(String, nullable=True, index=True)
    phylum_id = Column(String, nullable=True, index=True)
    class_id = Column(String, nullable=True, index=True)
    order_id = Column(String, nullable=True, index=True)
    family_id = Column(String, nullable=True, index=True)
    cofactors = Column(String, nullable=True, index=True)


class Reaction(Base):
    """Child of Entry (one or more reactions per entry)"""

    __tablename__ = "reaction"

    id = Column(Integer, primary_key=True)
    entry_id = Column(String, ForeignKey("entry.accession"))
    entry = relationship("Entry", back_populates="reactions")

    description = Column(Text, nullable=True)
    reaction_smarts = Column(Text)
    rhea_id = Column(String, nullable=True, index=True)
    ec_id = Column(String, nullable=True, index=True)

    example_reactions = relationship(
        "ExampleReaction", back_populates="reaction", cascade="all, delete-orphan"
    )


class ExampleReaction(Base):
    """Child of Reaction (one or more example reactions per reaction)"""

    __tablename__ = "example_reaction"

    id = Column(Integer, primary_key=True)
    reaction_id = Column(Integer, ForeignKey("reaction.id"))
    reaction = relationship("Reaction", back_populates="example_reactions")

    smiles_substrate = Column(Text)
    is_intermediate = Column(Boolean)

    products = relationship(
        "Product", back_populates="example_reaction", cascade="all, delete-orphan"
    )


class Product(Base):
    """Child of Example Reaction (one or more products per example reaction)"""

    __tablename__ = "product"

    id = Column(Integer, primary_key=True)
    example_reaction_id = Column(Integer, ForeignKey("example_reaction.id"))

    smiles_product = Column(Text)

    example_reaction = relationship("ExampleReaction", back_populates="products")
