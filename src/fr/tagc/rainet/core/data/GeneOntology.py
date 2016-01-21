
from sqlalchemy import Column, String
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base

from fr.tagc.rainet.core.data.ProteinGOAnnotation import ProteinGOAnnotation

##
# This class describe a gene ontology term

class GeneOntology( Base):
    __tablename__ = "GeneOntology"
    
    # The GO ID
    goID = Column( String, primary_key = True  )
    # The GO Term
    goTerm = Column( String )
    # The GO Aspect
    goAspect = Column( String)
    # Define the N-to-N relationship between GeneOntology and Protein
    annotatedProteins = relationship('Protein', secondary=ProteinGOAnnotation.__table__, backref="goAnnotations")
    
    def __init__(self, id, term, aspect):
        
        self.goID = id
        self.goTerm = term
        self.goAspect = aspect
        
    ##
    # Add an annotated protein to the list
    def add_annotated_protein(self, protein):
        
        if protein != None:
            self.annotatedProteins.append( protein)
    