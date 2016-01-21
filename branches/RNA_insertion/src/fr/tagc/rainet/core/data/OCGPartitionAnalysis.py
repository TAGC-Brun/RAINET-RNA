
from sqlalchemy import Column, String
from sqlalchemy.sql.schema import ForeignKey

from fr.tagc.rainet.core.data.PartitionAnalysis import PartitionAnalysis

# #
# This class describe a Paritition Analysis performed by OCG
#
class OCGPartitionAnalysis ( PartitionAnalysis ):
    __tablename__ = "OCGPartitionAnalysis"
    
    # The Partition Analysis name
    analysisName = Column( String, ForeignKey('PartitionAnalysis.analysisName'), primary_key=True)
    
    __mapper_args__ = {
        'polymorphic_identity':'OCGPartitionAnalysis',
    }
    
    def __init__(self, name, details):
        
        self.analysisName = name
        self.analysisDetails = details
    
    # #
    # Add a NetworkModule to the networkModules list
    #
    # @param module : NetworkModule - A NetworkModuleA ProteinIsoform instance
    def add_network_module( self, module ):
        
        if module != None:
            self.networkModules.append( module )
