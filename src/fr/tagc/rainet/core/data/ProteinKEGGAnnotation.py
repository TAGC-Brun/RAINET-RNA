
from sqlalchemy import Column, String, ForeignKey
from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from fr.tagc.rainet.core.util.exception.RainetException import RainetException

from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.Protein import Protein

# #
# This class describes a annotation of protein on KEGG Pathway
#
class ProteinKEGGAnnotation( Base ):
    __tablename__ = 'ProteinKEGGAnnotation'
    
    # The KEGG pathway annotation
    keggPathway_id = Column( String, ForeignKey( 'KEGGPathway.keggID', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)
    # The annotated protein
    protein_id = Column( String, ForeignKey( 'Protein.uniprotAC', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)

    # The cross reference
    crossReferenceID = Column( String )
    # Define the Composite PrimaryKey
#     __table_args__ = (
#         PrimaryKeyConstraint('protein_id', 'keggPathway_id', ),
#     )
    
    
    #
    # The constructor of the class
    #
    # @param cross_reference_id : string - The ID of the protein in BEGG DB notation
    # @param kegg_pathway_id : string - The ID of the KEGG pathway associated to the protein
    def __init__(self, cross_reference_id, kegg_pathway_id):
         
        self.crossReferenceID = cross_reference_id
         
        sql_session = SQLManager.get_instance().get_session()
         
        #=======================================================================
        # Search for the KEGGPathway corresponding to the given Kegg pathway ID
        #=======================================================================
        # -- make the query
        from fr.tagc.rainet.core.data.KEGGPathway import KEGGPathway
        kegg_pathway_list = sql_session.query( KEGGPathway).filter( KEGGPathway.keggID == kegg_pathway_id).all()
         
        # --Check if a single cross reference is found. If not, raise an issue
        kegg_pathway = None
        if kegg_pathway_list != None and len( kegg_pathway_list) > 0:
            if len( kegg_pathway_list) == 1:
                kegg_pathway = kegg_pathway_list[0]
            else:
                raise RainetException( "ProteinKEGGAnnotation.init : Abnormal number of KEGGPathway found for keggID = " + kegg_pathway_id + " : " + str( len( kegg_pathway_list))) 
        else:
            raise NotRequiredInstantiationException( "ProteinKEGGAnnotation.init : No KEGGPathway found for keggID = " + kegg_pathway_id)
 
        # -- Check if the KEGGPathway found is not None
        if kegg_pathway == None:
            raise RainetException( "ProteinKEGGAnnotation.init : returned cross reference is None for " + kegg_pathway_id)
         
        #=======================================================================
        # Get the ProteinCrossReference corresponding to the provided KEGG protein id
        # (cross_reference) in order to convert it to Uniprot AC
        #=======================================================================
        # -- make the query
        protein_crossref_list = sql_session.query( ProteinCrossReference).filter( ProteinCrossReference.sourceDB == "KEGG", 
                                                          ProteinCrossReference.crossReferenceID == cross_reference_id).all()
         
        # --Check if a single cross reference is found. If not, raise an issue
        if protein_crossref_list == None or len( protein_crossref_list) == 0:
            raise NotRequiredInstantiationException( "ProteinKEGGAnnotation.init : No ProteinCrossReference found for cross reference = " + cross_reference_id)

     
        #=======================================================================
        # Search from the Proteins corresponding to the cross references found
        # in the previous step. Then build the association between the Protein
        # and the KEGGPathway 
        #=======================================================================
        for cross_ref in protein_crossref_list:
            # -- make the query
            protein_list = sql_session.query( Protein).filter( Protein.uniprotAC == cross_ref.protein_id).all()
             
            # --Check if a single Protein is found. If not, raise an issue
            protein = None
            if protein_list != None and len( protein_list) > 0:
                if len( protein_list) == 1:
                    protein = protein_list[0]
                else:
                    raise RainetException( "ProteinKEGGAnnotation.init : Abnormal number of Protein found for cross reference = " + cross_ref.protein_id + " : " + str( len( protein_list))) 
            else:
                raise NotRequiredInstantiationException( "ProteinKEGGAnnotation.init : No Protein found for uniprotAC = " + cross_ref.protein_id)
     
            # -- Check if the Protein found is not None
            if protein == None:
                raise RainetException( "ProteinKEGGAnnotation.init : returned Protein is None for UniprotAC" + cross_ref.protein_id)
                 
            # -- Build the relation between the KEGGPathway and the Protein
            kegg_pathway.add_annotated_protein( protein)
            sql_session.add( kegg_pathway)
            sql_session.add( protein)
            
            # -- Search for the inserted annotation to add supplementary info
            annotation_list = sql_session.query( ProteinKEGGAnnotation).filter( ProteinKEGGAnnotation.keggPathway_id == kegg_pathway.keggID, ProteinKEGGAnnotation.protein_id == protein.uniprotAC).all()
             
            # --Check if a single ProteinKEGGAnnotation is found. If not, raise an issue
            # --If a single annotation is found, add the complementary info
            annotation = None
            if annotation_list != None and len( annotation_list) > 0:
                if len( annotation_list) == 1:
                    annotation = annotation_list[0]
                    annotation.crossReferenceID = cross_reference_id
                    sql_session.add( annotation)
                else:
                    raise RainetException( "ProteinKEGGAnnotation.init : Abnormal number of ProteinKEGGAnnotation found for uniprotAC= " + protein.uniprotAC + " and keggID=" + kegg_pathway.keggID + " : " + str( len( annotation_list))) 
            else:
                raise RainetException( "ProteinKEGGAnnotation.init : No ProteinKEGGAnnotation found for uniprotAC = " + protein.uniprotAC + " and keggID=" + kegg_pathway.keggID)

        # Raise the Exception to indicate the instance must not be inserted since it is automatically created by the KEGGPatway to Protein association done just before
        raise NotRequiredInstantiationException( "ProteinKEGGAnnotation.init : ProteinKEGGAnnotation objects do not have to be inserted by __init__ since they are created by KEGGPatway to Protein association table.")
        
            
    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
    
