
from sqlalchemy.orm import relationship
from sqlalchemy import Column, String

from fr.tagc.rainet.core.util.sql.Base import Base
from sqlalchemy.sql.schema import ForeignKey

# #
# This class describe a Partition Analysis performed by OCG
#
class PartitionAnalysis ( Base ):
    __tablename__ = "PartitionAnalysis"
    
    # The Partition Analysis name
    analysisName = Column( String, primary_key = True )
    # The partition details
    analysisDetails = Column( String )
    # The PPINetwork associated to this analysis
    networkName = Column( String, ForeignKey('PPINetwork.name'))
    # The list of network modules
    networkModules = relationship( 'NetworkModule')
    # Type of PartitionAnalysis
    type = Column( String )
    
    __mapper_args__ = {
        'polymorphic_identity':'PartitionAnalysis',
        'polymorphic_on':type
    }
    
    def __init__( self, name, details ):
        
        self.analysisName = name
        self.analysisDetails = details
    
    # #
    # Add a NetworkModule to the networkModules list
    #
    # @param module : NetworkModule - A NetworkModuleA ProteinIsoform instance
    def add_network_module( self, module ):
        
        if module != None:
            self.networkModules.append( module )
