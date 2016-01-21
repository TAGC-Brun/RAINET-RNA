
from sqlalchemy import Column, String, ForeignKey
from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from fr.tagc.rainet.core.util.exception.RainetException import RainetException

from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.util.log.Logger import Logger

# #
# This class describes a annotation of protein on KEGG Pathway
#
class ProteinReactomeAnnotation( Base ):
    __tablename__ = 'ProteinReactomeAnnotation'
    
    # The KEGG pathway annotation
    reatomePathway_id = Column( String, ForeignKey( 'ReactomePathway.reactomeID', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)
    # The annotated protein
    protein_id = Column( String, ForeignKey( 'Protein.uniprotAC', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)


    ##
    # The constructor of the class
    #
    # @param uniprot_ac : string - The ID of the protein in uniprot
    # @param reactome_pathway_id : string - The ID of the reactome pathway associated to the protein
    def __init__(self, protein_id, reactome_pathway_id):
         
        sql_session = SQLManager.get_instance().get_session()
         
        #=======================================================================
        # Search for the Reactome Pathway corresponding to the given Reactome pathway ID
        #=======================================================================
        # -- make the query
        from fr.tagc.rainet.core.data.ReactomePathway import ReactomePathway
        reactome_pathway_list = sql_session.query( ReactomePathway).filter( ReactomePathway.reactomeID == reactome_pathway_id).all()
         
        # --Check if a single cross reference is found. If not, raise an issue
        reactome_pathway = None
        if reactome_pathway_list != None and len( reactome_pathway_list) > 0:
            if len( reactome_pathway_list) == 1:
                reactome_pathway = reactome_pathway_list[0]
            else:
                raise RainetException( "ProteinReactomeAnnotation.init : Abnormal number of ReactomePathway found for reactomeID = " + reactome_pathway_id + " : " + str( len( reactome_pathway_list))) 
        else:
            raise NotRequiredInstantiationException( "ProteinReactomeAnnotation.init : No ReactomePathway found for reactomeID = " + reactome_pathway_id)
 
        # -- Check if the ReactomePathway found is not None
        if reactome_pathway == None:
            raise RainetException( "ProteinKEGGAnnotation.init : returned cross reference is None for " + reactome_pathway_id)

     
        #=======================================================================
        # Search from the Proteins corresponding to the uniprotAC
        #=======================================================================
        protein_list = sql_session.query( Protein).filter( Protein.uniprotAC == protein_id).all()
             
        # --Check if a single Protein is found. If not, raise an issue
        protein = None
        if protein_list != None and len( protein_list) > 0:
            if len( protein_list) == 1:
                protein = protein_list[0]
            else:
                raise RainetException( "ProteinReactomeAnnotation.init : Abnormal number of Protein found for uniprotAC= " + protein_id + " : " + str( len( protein_list))) 
        else:
            raise NotRequiredInstantiationException( "ProteinReactomeAnnotation.init : No Protein found for uniprotAC = " + protein_id)
 
        # -- Check if the Protein found is not None
        if protein == None:
            raise RainetException( "ProteinReactomeAnnotation.init : returned Protein is None for UniprotAC = " + protein_id)
             
        # -- Build the relation between the KEGGPathway and the Protein
        Logger.get_instance().debug(" Found protein = " + protein.uniprotAC)
        Logger.get_instance().debug(" Found pathway = " + reactome_pathway.reactomeID)
        reactome_pathway.add_annotated_protein( protein)
        sql_session.add( reactome_pathway)
        sql_session.add( protein)
        
        # Raise the Exception to indicate the instance must not be inserted since it is automatically created by the ReactomePatway to Protein association done just before
        raise NotRequiredInstantiationException( "ProteinReactomeAnnotation.init : ProteinReactomeAnnotation objects do not have to be inserted by __init__ since they are created by ReactomePatway to Protein association table.")
            
    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
    
