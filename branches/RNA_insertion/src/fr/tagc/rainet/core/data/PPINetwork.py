
from sqlalchemy import Column, String
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

# #
# This class describe a Protein-to-Protein Interaction Network
#
# this class has a one-to-many relationship with the class PPINetworkInteraction
#
class PPINetwork( Base ):
    __tablename__ = 'PPINetwork'

    # The network name
    name = Column( String, primary_key = True )
    # The list of ParitionAnalysis of the PPINetwork
    partitionAnalysis = relationship( 'PartitionAnalysis', backref="network")
    # The list of interaction in the PPI network
#     interactions = relationship( 'PPINetworkInteraction', backref = 'network' )

    
    # #
    # The Protein constructor
    # 
    # @param name : string - The name of the network
    def __init__( self, name):
        
        self.name = name
        
    # #
    # Add the object to SQLAlchemy session
    def add_to_session( self ):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self )

    # #
    # Add a PartitionAnalysis to the PPINetwork
    def add_partition_analysis(self, partition_analysis):
        self.partitionAnalysis.append( partition_analysis)
        
#     # #
#     # Add an interaction to the network
#     #
#     # @param interaction : string - A PPINetworkInteraction instance
#     def add_interaction( self, interaction ):
#         
#         if interaction != None:
#             self.interactions.append( interaction )
#             print "Adding interaction to PPINetwork"
        