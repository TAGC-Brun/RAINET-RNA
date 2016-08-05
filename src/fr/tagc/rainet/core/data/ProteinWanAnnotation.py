
from sqlalchemy import Column, String, ForeignKey
from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from fr.tagc.rainet.core.util.exception.RainetException import RainetException

from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.SynonymGeneSymbol import SynonymGeneSymbol

# #
# This class describes a annotation of protein on Wan clusters
#
class ProteinWanAnnotation( Base ):
    __tablename__ = 'ProteinWanAnnotation'
    
    # The Wan cluster annotation
    wan_cluster_id = Column( String, ForeignKey( 'WanCluster.clusterID', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)
    # The annotated protein
    protein_id = Column( String, ForeignKey( 'Protein.uniprotAC', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)  
    
    #
    # The constructor of the class
    #
    # @param wan_cluster_id : string - The ID of the Wan cluster associated to the protein
    def __init__(self, wan_cluster_id, protein_id):
                
        sql_session = SQLManager.get_instance().get_session()
         
        #=======================================================================
        # Search for the Wan cluster corresponding to the given Wan cluster ID
        #=======================================================================
        # -- make the query
        from fr.tagc.rainet.core.data.WanCluster import WanCluster
        clusters_list = sql_session.query( WanCluster).filter( WanCluster.clusterID == wan_cluster_id).all()
         
        # --Check if a single cross reference is found. If not, raise an issue
        cluster_id = None
        if clusters_list != None and len( clusters_list) > 0:
            if len( clusters_list) == 1:
                cluster_id = clusters_list[0]
            else:
                raise RainetException( "ProteinWanAnnotation.init : Abnormal number of Wan clusters found for cluster_id = " + cluster_id + " : " + str( len( clusters_list))) 
        else:
            raise NotRequiredInstantiationException( "ProteinWanAnnotation.init : No Wan cluster found for cluster id = " + cluster_id)
 
        if cluster_id == None:
            raise RainetException( "ProteinWanAnnotation.init : returned cross reference is None for " + cluster_id)
         
        #=======================================================================
        # Search the protein object. Then build the association between the Protein and the WanCluster 
        #=======================================================================
        # -- make the query
        protein_list = sql_session.query( Protein).filter( Protein.uniprotAC == protein_id).all()
         
        # --Check if a single Protein is found. If not, raise an issue
        protein = None
        if protein_list != None and len( protein_list) > 0:
            if len( protein_list) == 1:
                protein = protein_list[0]
            else:
                raise RainetException( "ProteinWanAnnotation.init : Abnormal number of Protein found for = " +protein_id + " : " + str( len( protein_list))) 
        else:
            raise NotRequiredInstantiationException( "ProteinWanAnnotation.init : No Protein found for uniprotAC = " + protein_id)
 
        # -- Check if the Protein found is not None
        if protein == None:
            raise RainetException( "ProteinWanAnnotation.init : returned Protein is None for UniprotAC" + protein_id)
             
        # -- Build the relation between the Wan cluster and the Protein
        cluster_id.add_annotated_protein( protein)
        sql_session.add( cluster_id)
        sql_session.add( protein)
            
        # Raise the Exception to indicate the instance must not be inserted since it is automatically created 
        raise NotRequiredInstantiationException( "ProteinWanAnnotation.init : ProteinWanAnnotation objects do not have to be inserted by __init__ since they are created by WanCluster to Protein association table.")
        
            
    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
    
