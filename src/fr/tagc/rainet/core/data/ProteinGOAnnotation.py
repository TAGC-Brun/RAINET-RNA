
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
class ProteinGOAnnotation( Base ):
    __tablename__ = 'ProteinGOAnnotation'
    
    # The GO term
    geneOntology_id = Column( String, ForeignKey( 'GeneOntology.goID', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)
    # The annotated protein
    protein_id = Column( String, ForeignKey( 'Protein.uniprotAC', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)
    # The annotation evidence code
    evidenceCode = Column( String)
    
    ##
    # The constructor of the class
    #
    # @param uniprot_ac : string - The ID of the protein in uniprot
    # @param gene_ontology_id : string - The ID of the GeneOntology associated to the protein
    def __init__(self, source_sb, protein_id, gene_ontology_id, evidence_code):
         
        sql_session = SQLManager.get_instance().get_session()
         
        #=======================================================================
        # Search for the GeneOntology corresponding to the given GeneOntology ID
        #=======================================================================
        # -- make the query
        from fr.tagc.rainet.core.data.GeneOntology import GeneOntology
        gene_ontology_list = sql_session.query( GeneOntology).filter( GeneOntology.goID == gene_ontology_id).all()
         
        # --Check if a single cross reference is found. If not, raise an issue
        gene_ontology = None
        if gene_ontology_list != None and len( gene_ontology_list) > 0:
            if len( gene_ontology_list) == 1:
                gene_ontology = gene_ontology_list[0]
            else:
                raise RainetException( "ProteinGOAnnotation.init : Abnormal number of GeneOntology found for goID = " + gene_ontology_id + " : " + str( len( gene_ontology_list))) 
        else:
            raise NotRequiredInstantiationException( "ProteinGOAnnotation.init : No GeneOntology found for goID = " + gene_ontology_id)
 
        # -- Check if the GeneOntology found is not None
        if gene_ontology == None:
            raise RainetException( "ProteinGOAnnotation.init : returned cross reference is None for " + gene_ontology_id)

     
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
                raise RainetException( "ProteinGOAnnotation.init : Abnormal number of Protein found for uniprotAC= " + protein_id + " : " + str( len( protein_list))) 
        else:
            raise NotRequiredInstantiationException( "ProteinGOAnnotation.init : No Protein found for uniprotAC = " + protein_id)
 
        # -- Check if the Protein found is not None
        if protein == None:
            raise RainetException( "ProteinGOAnnotation.init : returned Protein is None for UniprotAC = " + protein_id)
             
        # -- Build the relation between the KEGGPathway and the Protein
        gene_ontology.add_annotated_protein( protein)
        sql_session.add( gene_ontology)
        sql_session.add( protein)
        
        #=======================================================================
        # Search for the inserted annotation to add supplementary info
        #=======================================================================
        annotation_list = sql_session.query( ProteinGOAnnotation).filter( ProteinGOAnnotation.geneOntology_id == gene_ontology.goID, ProteinGOAnnotation.protein_id == protein.uniprotAC).all()
         
        # --Check if a single ProteinGOAnnotation is found. If not, raise an issue
        # --If a single annotation is found, add the complementary info
        annotation = None
        if annotation_list != None and len( annotation_list) > 0:
            if len( annotation_list) == 1:
                annotation = annotation_list[0]
                annotation.evidenceCode = evidence_code
                sql_session.add( annotation)
            else:
                raise RainetException( "ProteinGOAnnotation.init : Abnormal number of ProteinGOAnnotation found for uniprotAC= " + protein_id + " and go=" + gene_ontology_id + " : " + str( len( annotation_list))) 
        else:
            raise RainetException( "ProteinGOAnnotation.init : No ProteinGOAnnotation found for uniprotAC = " + protein_id + " and go=" + gene_ontology_id)
        
        # Raise the Exception to indicate the instance must not be inserted since it is automatically created by the ReactomePatway to Protein association done just before
        raise NotRequiredInstantiationException( "ProteinGOAnnotation.init : ProteinGOAnnotation objects do not have to be inserted by __init__ since they are created by GeneOntology to Protein association table.")
            
    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
    
