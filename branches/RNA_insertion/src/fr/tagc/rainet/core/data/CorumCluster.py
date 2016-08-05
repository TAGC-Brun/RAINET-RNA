
from sqlalchemy import Column, String
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.ProteinCorumAnnotation import ProteinCorumAnnotation

# #
# This class describes a Corum cluster
#
class CorumCluster( Base ):
    __tablename__ = 'CorumCluster'
    
    # The pathway ID
    clusterID = Column( String, primary_key = True  )
    # Cluster description
    clusterName = Column( String)
    # Method used for identification
    clusterMethod = Column( String) 

    # Define the N-to-N relationship between CorumAnnotation and Protein
    annotatedProteins = relationship('Protein', secondary=ProteinCorumAnnotation.__table__, backref="CorumAnnotations")
    
    #
    # The constructor of the class
    #
    # @param cluster_id : string - The ID of the Corum cluster
    def __init__(self, cluster_id, cluster_name, cluster_method):
        
        self.clusterID = cluster_id
        self.clusterName = cluster_name
        self.clusterMethod = cluster_method
            
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
