
from sqlalchemy import Column, String
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.ProteinReactomeAnnotation import ProteinReactomeAnnotation

# #
# This class describes a Reactome pathway
#
class ReactomePathway( Base ):
    __tablename__ = 'ReactomePathway'
    
    # The reactome pathway ID
    reactomeID = Column( String, primary_key = True  )
    # The reactome pathway name
    reactomeName = Column( String )
    # Define the N-to-N relationship between Reactome and Protein
    annotatedProteins = relationship('Protein', secondary=ProteinReactomeAnnotation.__table__, backref="reactomeAnnotations")
    
    #
    # The constructor of the class
    #
    # @param reactome_id : string - The ID of the reactome pathway
    # @param reactome_name : string - The name of the reactome pathway
    def __init__(self, reactome_id, reactome_name):
        
        self.reactomeID = reactome_id
        self.reactomeName = str( reactome_name)
    
    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        identical_reactome_pathway = sql_session.query( ReactomePathway).filter( ReactomePathway.reactomeID == self.reactomeID).all()
        if identical_reactome_pathway == None or len( identical_reactome_pathway) == 0:
            sql_session.add( self)
    
    ##
    # Add an annotated protein to the list
    def add_annotated_protein(self, protein):
        
        if protein != None:
            self.annotatedProteins.append( protein)
