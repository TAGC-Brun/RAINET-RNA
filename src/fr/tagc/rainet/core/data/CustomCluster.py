
from sqlalchemy import Column, String
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.ProteinCustomAnnotation import ProteinCustomAnnotation

# #
# This class describes a Custom cluster
#
class CustomCluster( Base ):
    __tablename__ = 'CustomCluster'
    
    # The pathway ID
    clusterID = Column( String, primary_key = True  )
    # Define the N-to-N relationship between CustomAnnotation and Protein
    annotatedProteins = relationship('Protein', secondary=ProteinCustomAnnotation.__table__, backref="CustomAnnotations")
    
    #
    # The constructor of the class
    #
    # @param cluster_id : string - The ID of the Custom cluster
    def __init__(self, cluster_id):
        
        self.clusterID = cluster_id
            
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
