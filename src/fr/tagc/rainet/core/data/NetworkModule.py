from sqlalchemy import Column, String
from sqlalchemy.orm import relationship, backref
from sqlalchemy.sql.schema import ForeignKey

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.data.ProteinNetworkModule import ProteinNetworkModule


##
# This class describes a Network Module (a group of interacting proteins)
# as produced by a network partition analysis like OCG
class NetworkModule ( Base ):
    __tablename__ = "NetworkModule"
    
    # The module ID
    moduleID = Column( String, primary_key = True )
    # The origin PartitionAnalysis
    partitionAnalysis_id = Column( String, ForeignKey('PartitionAnalysis.analysisName'))
    # Define the N-to-N relationship between NetworkModule and Protein
    associatedProteins = relationship('Protein', secondary=ProteinNetworkModule.__table__, backref="networkModules")
    # The list of module annotations
    associatedAnnotations = relationship( 'NetworkModuleAnnotation', backref= "networkModule")
    
    #
    ## The constructor
    def __init__(self, module_id):
        
        self.moduleID = module_id
    
        
    ##
    # Add an associated protein to the list
    def add_associated_protein(self, protein):
        
        if protein != None:
            self.associatedProteins.append( protein)

    ##
    # Add an associated annotation to the list
    def add_associated_annotation(self, association):
        
        if association != None:
            self.associatedAnnotations.append( association)