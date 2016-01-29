
import re

from sqlalchemy import or_,and_
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger

from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.pattern.PatternUtil import PatternUtil
from fr.tagc.rainet.core.util.sql.SQLUtil import SQLUtil

from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference

from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from sqlalchemy.orm.base import instance_dict


# #
# This class describes a RNA models obtained from querying Ensembl BioMart
#
class ProteinRNAInteraction( Base ):
    __tablename__ = 'ProteinRNAInteraction'

    # The RNA transcript ID 
    transcriptID = Column( String, ForeignKey( 'RNA.transcriptID'), primary_key=True)
    # The Protein transcript ID
    proteinID = Column( String, ForeignKey( 'Protein.uniprotAC'), primary_key=True)

    interactionScore = Column( String )
    
    proteins = relationship( "Protein" )


    # #
    # The ProteinRNAInteraction constructor.
    #
    def __init__( self, interactors, interaction_score):

        sql_session = SQLManager.get_instance().get_session()

        #=======================================================================
        # Fill the main RNA variables
        #=======================================================================

        try:
            self.interactionScore = float(interaction_score)
        except ValueError as ve:
            raise RainetException( "ProteinRNAInteraction.__init__ : The value of interaction score is not a float: " + str( interaction_score ), ve )

        #=======================================================================
        # Parse interactors
        # dr note: in the initial catRAPID input, both interactors are in the same TSV column. E.g.  ENSP00000320917_ENST00000389594
        # Protein is always on left side, RNA in the right side. IDs should not contain "_".
        #=======================================================================

        spl = interactors.split("_")
        if len(spl) == 2:
            protein_id = spl[0]
            transcript_id = spl[1]
        else:
            raise RainetException( "ProteinRNAInteraction.__init__ : The interactor string could not be parsed: " + str( interactors ))
            
        #=======================================================================
        # Query Protein object using cross references
        #=======================================================================
 
        from fr.tagc.rainet.core.data.Protein import Protein
        Logger.get_instance().debug( "\nSearching proteins" )

        #From Cross ref ID, search primary ID using cross reference table
        uniprot_ac = sql_session.query( ProteinCrossReference.protein_id ).filter( and_(ProteinCrossReference.sourceDB == "Ensembl_PRO", ProteinCrossReference.crossReferenceID == protein_id ) ).first()
        uniprot_ac = str(uniprot_ac[0])
        # dr: this needs to be properly coded, getting exceptions, using constants etc. See Lionel's find_protein

        #Get Protein instance using primary ID
        protein_list = sql_session.query( Protein ).filter( or_( Protein.uniprotAC == uniprot_ac)).all()
        if protein_list != None and len( protein_list) == 1 :
            protein = protein_list[0]
            self.proteinID = protein
        else:
            raise RainetException( "ProteinRNAInteraction.__init__ : Could not retrieve protein object: " + str( protein_id ))
        ###throw warning if protein not found, instead of exception
        
#        self.proteinID = self.find_protein( protein_id )

 
        #=======================================================================
        # Query RNA object
        #=======================================================================

        # Retrieve the list of RNAs corresponding to the provided accession number
        from fr.tagc.rainet.core.data.RNA import RNA

        RNA_list = sql_session.query( RNA ).filter( or_( RNA.transcriptID == transcript_id)).all()

        if RNA_list != None and len( RNA_list) == 1 :
            RNA = RNA_list[0]
            self.transcriptID = RNA
        else:
            raise RainetException( "ProteinRNAInteraction.__init__ : Could not retrieve RNA object: " + str( transcript_id ))
        ###throw warning if RNA not found, instead of exception

         
         
        self.add_to_session()


    # #
    # Search for a protein ID (uniprotAC or other) in the ids of the provided id_list
    # 
    def find_protein( self, protein_id ):
        
        #=======================================================================
        # Search for a UniprotAC in the provided IDs
        #=======================================================================
    
        Logger.get_instance().debug( "|--Find_proteins :Trying crossreferences ---" )
        # Initialize uniprot_ac_list
        uniprot_ac_list = []
        
        # Get a SQLalchemy session
        sql_session = SQLManager.get_instance().get_session()

        print (protein_id)

        uniprot_ac = sql_session.query( ProteinCrossReference.protein_id ).filter( and_(ProteinCrossReference.sourceDB == "Ensembl_PRO", ProteinCrossReference.crossReferenceID == protein_id ) ).all()
        return (str(uniprot_ac[0])) #should I add / return object or ID?
    
#         # Parse the cross reference regular expression to find matching ID that will be
#         # converted to uniprotAC thanks to cross references present in database
#         for crossref_source, crossref_pattern in DataConstants.PRI_ID_CROSSREF_REGEX_DICT.items():
#             Logger.get_instance().debug( "|--|--Cross reference = " + crossref_source )
#             # Compile the regular expression
#             crossref_compiled_pattern = re.compile( crossref_pattern )
#             # Search matching string (if any) in the provided strings 
#             matched_crossref_list = PatternUtil.find_groups_in_list( crossref_compiled_pattern, id_list )
#             Logger.get_instance().debug( "|--|--|--Cross reference list found = " + str( matched_crossref_list ) )
#             # If matching string exists, look in database if corresponding crossreference is found
#             # If so, add the equivalent uniprotAC to the final list.
#             if matched_crossref_list != None and len( matched_crossref_list ) > 0:
#                 for matched_crossref in matched_crossref_list:
#                     db_regex = DataConstants.INTERACTOME_ID_CROSSREF_DB_REGEX_DICT[ crossref_source].replace( "VALUE", matched_crossref )
#                     uniprot_ac = sql_session.query( ProteinCrossReference.protein_id ).filter( ProteinCrossReference.sourceDB == crossref_source, ProteinCrossReference.crossReferenceID.like( db_regex ) ).all()
#                     if uniprot_ac != None and len( uniprot_ac) >= 1:
#                         uniprot_ac_list.extend( uniprot_ac[0] )
#         
#         Logger.get_instance().debug( "|--Found unitproAC= " + str( uniprot_ac_list ) )
#         
#         #=======================================================================
#         # Return the protein with the found uniprot_ac
#         #=======================================================================
#         try:
#             return SQLUtil.find_unique_protein_from_uniprotAC_list( uniprot_ac_list )
#         except RainetException as e:
#             Logger.get_instance().debug( e.to_string() ) 
#             return None

    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
