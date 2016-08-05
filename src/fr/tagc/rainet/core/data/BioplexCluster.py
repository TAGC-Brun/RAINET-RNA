
from sqlalchemy import Column, String
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.ProteinBioplexAnnotation import ProteinBioplexAnnotation

# #
# This class describes a Bioplex cluster
#
class BioplexCluster( Base ):
    __tablename__ = 'BioplexCluster'
    
    # The pathway ID
    bioplexID = Column( String, primary_key = True  )
    # Define the N-to-N relationship between Bioplex cluster and Protein
    annotatedProteins = relationship('Protein', secondary=ProteinBioplexAnnotation.__table__, backref="bioplexAnnotations")
    
    #
    # The constructor of the class
    #
    # @param cluster_id : string - The ID of the Bioplex cluster
    def __init__(self, cluster_id):
        
        self.bioplexID = cluster_id
            
    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
    
    ##
    # Add an annotated protein to the list
    def add_annotated_protein(self, protein):
        
        if protein != None:
            self.annotatedProteins.append( protein)
