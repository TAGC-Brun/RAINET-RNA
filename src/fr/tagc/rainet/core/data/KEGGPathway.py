
from sqlalchemy import Column, String
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.ProteinKEGGAnnotation import ProteinKEGGAnnotation

# #
# This class describes a KEGG Pathway
#
class KEGGPathway( Base ):
    __tablename__ = 'KEGGPathway'
    
    # The pathway ID
    keggID = Column( String, primary_key = True  )
    # The pathway name
    keggName = Column( String )
    # Define the N-to-N relationship between KEGGPathway and Protein
    annotatedProteins = relationship('Protein', secondary=ProteinKEGGAnnotation.__table__, backref="keggAnnotations")
    
    #
    # The constructor of the class
    #
    # @param kegg_id : string - The ID of the KEGG Pathway
    # @param kegg_name : string - The name of the pathway
    def __init__(self, kegg_id, kegg_name):
        
        #  
        self.keggID = kegg_id
        
        # The real pathway name is the first part of the name (in the file the
        # organism is added like :"Glycolysis / Gluconeogenesis - Homo sapiens (human)"
        if kegg_name != None:
            token_list = kegg_name.split("-")
            if token_list != None and len(token_list) >= 1:
                self.keggName = token_list[0]
    
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
