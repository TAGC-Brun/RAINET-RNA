
from sqlalchemy import Column, String, ForeignKey, PrimaryKeyConstraint, ForeignKeyConstraint
from fr.tagc.rainet.core.util.sql.Base import Base

##
# This class describe a gene symbol that is a synonym of a main gene symbol of a protein
# this class has a many-to-one relationship with the class Protein
#
class SynonymGeneSymbol( Base):
    __tablename__ = 'SynonymGeneSymbol'
    
    # The protein associated to the related GeneSymbol
    protein_id = Column( String)
    # The related GeneSymbol 
    uniprotGeneSymbol_id = Column( String)
    # The synonym symbol
    uniprotSynonymGeneSymbol = Column( String)
    # Define the Composite PrimaryKey and ForeignKeys
    __table_args__ = (
        PrimaryKeyConstraint('protein_id', 'uniprotGeneSymbol_id', 'uniprotSynonymGeneSymbol'),
        ForeignKeyConstraint( ['protein_id', 'uniprotGeneSymbol_id'],
                            ['GeneSymbol.protein_id', 'GeneSymbol.uniprotGeneSymbol'])
                            ,)
