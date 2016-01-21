
from sqlalchemy import Column, String, Integer, ForeignKey,PrimaryKeyConstraint
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base

##
# This class describe a gene symbol associated to a protein
# this class has a many-to-one relationship with the class Protein
#
class GeneSymbol( Base):
    __tablename__ = 'GeneSymbol'
    
    # The base Protein
    protein_id = Column(String, ForeignKey('Protein.uniprotAC'))
    # The synonym symbol
    uniprotGeneSymbol = Column( String)
    # Define the Composite PrimaryKey
    __table_args__ = (
        PrimaryKeyConstraint('protein_id', 'uniprotGeneSymbol'),
    )
    # The list of Uniprot SynonymGeneSymbol
    uniprotSynonymGeneSymbols = relationship( 'SynonymGeneSymbol', backref = 'geneSymbol' )
    
    ##
    # Add a SynonymGeneSymbol to the list
    def add_synonym(self, synonym):
        
        if synonym != None:
            self.uniprotSynonymGeneSymbols.append( synonym)