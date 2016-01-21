import re

from sqlalchemy import Column, String, Float, PrimaryKeyConstraint, ForeignKeyConstraint
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base

from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.pattern.PatternUtil import PatternUtil
from fr.tagc.rainet.core.util.sql.SQLUtil import SQLUtil
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference

# #
# This class describes a Protein to Protein interaction
#
class ProteinInteraction( Base ):
    __tablename__ = 'ProteinInteraction'
    
    # The first interactor
    interactorA_id = Column( String )
    interactorA = relationship( 'Protein', foreign_keys = interactorA_id )
    # The first interactor
    interactorB_id = Column( String )
    interactorB = relationship( 'Protein', foreign_keys = interactorB_id )
    # The method used to detect the interaction
    detectionMethod = Column( String )
    # The pubmedID of the publication reporting the interaction
    pubmedID = Column( String )
    # The type of the interaction
    interactionType = Column( String )
    # The database reporting where the interaction was found
    sourceDB = Column( String )
    # The confidence score associated to the interaction
    confidenceScore = Column( Float )
    # Define the Composite PrimaryKey and ForeignKeys
    __table_args__ = ( 
        PrimaryKeyConstraint( 'interactorA_id', 'interactorB_id', 'detectionMethod', 'pubmedID' ),
        ForeignKeyConstraint( ['interactorA_id', 'interactorB_id'],
                            ['Protein.uniprotAC', 'Protein.uniprotAC'] )
                            , )
    
    # #
    # The constructor
    #
    def __init__( self, id_A, id_B, alt_id_A, alt_id_B, method, pubmed_string, interaction_type, source_db, score_string ):
        
        Logger.get_instance().debug( "\nSearching proteins" )
        protein_A = self.find_protein( [id_A, alt_id_A] )
        Logger.get_instance().debug( "|--proteinA found = " + str( protein_A) )
        
        protein_B = self.find_protein( [id_B, alt_id_B] )
        Logger.get_instance().debug( "|--proteinB found = " + str( protein_B ) )

        # If both interacting proteins have been found
        if protein_A != None and protein_B != None:
            self.pubmedID = PatternUtil.find_single_group_in_string( DataConstants.INTERACTOME_PUBMED_REGEX, pubmed_string.lower(), False )
#             self.detectionMethod = PatternUtil.find_single_group_in_string( DataConstants.INTERACTOME_METHOD_REGEX, method, False )
#             self.interactionType = PatternUtil.find_single_group_in_string( DataConstants.INTERACTOME_TYPE_REGEX, interaction_type, False )
#             self.sourceDB = PatternUtil.find_single_group_in_string( DataConstants.INTERACTOME_SOURCEDB_REGEX, source_db, False )
            self.detectionMethod= method
            self.interactionType = interaction_type
            self.sourceDB = source_db
            self.confidenceScore = self.get_score( score_string )
            Logger.get_instance().debug( "--->score is " + str( self.confidenceScore))
            
            # Test if a similar interaction already exists (may be possible due to isoform of proteins)
            sql_session = SQLManager.get_instance().get_session()
            if self.search_similar_interaction( protein_A, protein_B, sql_session ) != None:
                raise NotRequiredInstantiationException( "ProteinInteraction: A similar interaction exists in database." )
            
            # Add interaction to interacting proteins
            self.interactorA = protein_A
            self.interactorB = protein_B
            sql_session.add( protein_A )
            sql_session.add( protein_B )
        else:
            raise NotRequiredInstantiationException( "ProteinInteraction: At least one interacting protein was not found in database." )

    # #
    # Search for a protein ID (uniprotAC or other) in the ids of the provided id_list
    # 
    def find_protein( self, id_list ):
        
        #=======================================================================
        # Search for a UniprotAC in the provided IDs
        #=======================================================================
        uniprotkb_pattern = re.compile( DataConstants.INTERACTOME_ID_UNIPROTKB_REGEX )
        uniprot_ac_list = PatternUtil.find_groups_in_list( uniprotkb_pattern, id_list )
        # -- If the uniprotAC contains a "-", it is a isoform uniprotAC so only the base
        # of the AC must be kept
        
        # Search for other ID type using also cross references if no uniprot AC was found before
        if uniprot_ac_list == None or len ( uniprot_ac_list ) == 0:
            
            Logger.get_instance().debug( "|--Find_proteins :Trying crossrefences ---" )
            # Initialize uniprot_ac_list
            uniprot_ac_list = []
            
            # Get a SQLalchemy session
            sql_session = SQLManager.get_instance().get_session()
            
            # Parse the cross reference regular expression to find matching ID that will be
            # converted to uniprotAC thanks to cross references present in database
            for crossref_source, crossref_pattern in DataConstants.INTERACTOME_ID_CROSSREF_REGEX_DICT.items():
                Logger.get_instance().debug( "|--|--Cross reference = " + crossref_source )
                # Compile the regular expression
                crossref_compiled_pattern = re.compile( crossref_pattern )
                # Search matching string (if any) in the provided strings 
                matched_crossref_list = PatternUtil.find_groups_in_list( crossref_compiled_pattern, id_list )
                Logger.get_instance().debug( "|--|--|--Cross reference list found = " + str( matched_crossref_list ) )
                # If matching string exists, look in database if corresponding crossreference is found
                # If so, add the equivalent uniprotAC to the final list.
                if matched_crossref_list != None and len( matched_crossref_list ) > 0:
                    for matched_crossref in matched_crossref_list:
                        db_regex = DataConstants.INTERACTOME_ID_CROSSREF_DB_REGEX_DICT[ crossref_source].replace( "VALUE", matched_crossref )
                        uniprot_ac = sql_session.query( ProteinCrossReference.protein_id ).filter( ProteinCrossReference.sourceDB == crossref_source, ProteinCrossReference.crossReferenceID.like( db_regex ) ).all()
                        if uniprot_ac != None and len( uniprot_ac) >= 1:
                            uniprot_ac_list.extend( uniprot_ac[0] )
        
        Logger.get_instance().debug( "|--Found unitproAC= " + str( uniprot_ac_list ) )
        
        #=======================================================================
        # Return the protein with the found uniprot_ac
        #=======================================================================
        try:
            return SQLUtil.find_unique_protein_from_uniprotAC_list( uniprot_ac_list )
        except RainetException as e:
            Logger.get_instance().debug( e.to_string() ) 
            return None
        

    # #
    # Retrieve the interaction score in the given string
    #
    # @param score_string : string - The string containing the score like "intact-miscore:0.96"
    #
    # @return float - The interaction score or -1 if an issue occurred
    #
    def get_score( self, score_string):
        
        compiled_pattern = re.compile( DataConstants.INTERACTOME_SCORE_REGEX )
        score = PatternUtil.find_groups_in_string( compiled_pattern, score_string )
        
        Logger.get_instance().debug( "|--|--|-score groups = " + str( score ) )
        if score == None:
            return -1
        elif len( score) >= 1:
            score_number = score[0]
        else:
            return -1
        
        try:
            return float( score_number )
        except ValueError:
            Logger.get_instance().warning( "ProteinInteraction.get_score : score is not a number : " + score_number )
            return -1    
        
    ##
    # Search for a ProteinInteraction between protein_a and protein_b with detectionMethod and PubmedID 
    # equals to those inserted in the current ProteinInteraction. The search is done for protein in
    # both directions.
    #
    # @param protein_a : Protein - The first protein
    # @param protein_b : Protein - The second protein
    # @param sql_session : SQLSession - The current SQLAlchemy session
    def search_similar_interaction( self, protein_a, protein_b, sql_session ):
        
        # Look at an interaction from A to B
        a_to_b_list = sql_session.query( ProteinInteraction ).filter( ProteinInteraction.interactorA == protein_a,
                                                       ProteinInteraction.interactorB == protein_b,
                                                       ProteinInteraction.detectionMethod == self.detectionMethod,
                                                       ProteinInteraction.pubmedID == self.pubmedID ).all()
                                                       
        # Check the result: if a single interaction is found, return it. If several are found, raise an exception                                
        if a_to_b_list != None and len( a_to_b_list ) > 0:
            if len( a_to_b_list ) == 1:
                return a_to_b_list[0]
            else:
                raise RainetException( "ProteinInteraction.__init__ : Found several protein interaction for interactor A = " + protein_a.uniprotAC + 
                                      " /interactor B = " + protein_b.uniprotAC + " /detection method = " + self.detectionMethod + " /pubmedID = " + self.pubmedID )
        
        # Look at an interaction from B to A
        b_to_a_list = sql_session.query( ProteinInteraction ).filter( ProteinInteraction.interactorA == protein_b,
                                                       ProteinInteraction.interactorB == protein_a,
                                                       ProteinInteraction.detectionMethod == self.detectionMethod,
                                                       ProteinInteraction.pubmedID == self.pubmedID ).all()
        
        # Check the result: if a single interaction is found, return it. If several are found, raise an exception
        if b_to_a_list != None and len( b_to_a_list ) > 0:
            if len( b_to_a_list ) == 1:
                return b_to_a_list[0]
            else:
                raise RainetException( "ProteinInteraction.__init__ : Found several protein interaction for interactor A = " + protein_b.uniprotAC + 
                                      " /interactor B = " + protein_a.uniprotAC + " /detection method = " + self.detectionMethod + " /pubmedID = " + self.pubmedID )
                
        return None
    
    # #
    # Add the object to SQLAlchemy session
    def add_to_session( self ):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self )
