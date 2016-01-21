from sqlalchemy import Column, String, PrimaryKeyConstraint, ForeignKey, ForeignKeyConstraint
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base

from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.sql.SQLUtil import SQLUtil
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException

from fr.tagc.rainet.core.data.PPINetwork import PPINetwork

#from fr.tagc.rainet.core.data.PPINetwork import PPINetwork

# #
# This class describes a Protein to Protein interaction in a PPI Network
#
class PPINetworkInteraction( Base ):
    __tablename__ = 'PPINetworkInteraction'
    
    # The first interactor
    interactorA_id = Column( String )
    interactorA = relationship( 'Protein', foreign_keys = interactorA_id )
    # The first interactor
    interactorB_id = Column( String )
    interactorB = relationship( 'Protein', foreign_keys = interactorB_id )
    # the name of the PPI network the interaction is part of
    network_id = Column( String)
    # Define the Composite PrimaryKey and ForeignKeys
    __table_args__ = ( 
        PrimaryKeyConstraint( 'interactorA_id', 'interactorB_id', 'network_id'),
        ForeignKeyConstraint( ['interactorA_id', 'interactorB_id'],
                            ['Protein.uniprotAC', 'Protein.uniprotAC'] )                      
                            , )
    
    
    # #
    # The constructor
    #
    def __init__( self, id_A, id_B, network_name="Default"):
        
        sql_session = SQLManager.get_instance().get_session()
        
        # Found PPINetwork the interaction is related to. If no network with the given
        # name exists, creates one
        network = sql_session.query( PPINetwork ).filter( PPINetwork.name == network_name).first()
        if network == None:
            network = PPINetwork( network_name)
        
        # Search for the proteins composing the interaction
        Logger.get_instance().debug( "\nSearching proteins" )
        protein_A = self.find_protein( [id_A] )
        Logger.get_instance().debug( "|--proteinA found = " + str( protein_A) )
        
        protein_B = self.find_protein( [id_B] )
        Logger.get_instance().debug( "|--proteinB found = " + str( protein_B ) )

        # If both interacting proteins have been found
        if protein_A != None and protein_B != None:
            
            # Test if a similar interaction already exists (may be possible due to isoform of proteins)
            if self.search_similar_interaction( protein_A, protein_B, sql_session ) != None:
                raise NotRequiredInstantiationException( "PPINetworkInteraction: A similar interaction exists in database." )
            
            # Add interaction to interacting proteins
            self.interactorA = protein_A
            self.interactorB = protein_B
            sql_session.add( protein_A )
            sql_session.add( protein_B )
#            # Add the interaction to the PPINetwork
#             network.add_interaction( self)
            self.network_id = network_name
            sql_session.add( network)
        else:
            raise NotRequiredInstantiationException( "PPINetworkInteraction: At least one interacting protein was not found in database." )

    # #
    # Search for a protein ID (uniprot) in the ids of the provided id_list
    # 
    def find_protein( self, uniprot_ac_list ):
        
        try:
            return SQLUtil.find_unique_protein_from_uniprotAC_list( uniprot_ac_list )
        except RainetException as e:
            Logger.get_instance().debug( e.to_string() ) 
            return None
        

    ##
    # Search for a PPINetworkInteraction between protein_a and protein_b 
    # equals to those inserted in the current PPINetworkInteraction. The search is done for protein in
    # both directions.
    #
    # @param protein_a : Protein - The first protein
    # @param protein_b : Protein - The second protein
    # @param sql_session : SQLSession - The current SQLAlchemy session
    def search_similar_interaction( self, protein_a, protein_b, sql_session ):
        
        # Look at an interaction from A to B
        a_to_b_list = sql_session.query( PPINetworkInteraction ).filter( PPINetworkInteraction.interactorA == protein_a,
                                                       PPINetworkInteraction.interactorB == protein_b).all()
                                                       
        # Check the result: if a single interaction is found, return it. If several are found, raise an exception                                
        if a_to_b_list != None and len( a_to_b_list ) > 0:
            if len( a_to_b_list ) == 1:
                return a_to_b_list[0]
            else:
                raise RainetException( "PPINetworkInteraction.__init__ : Found several protein interaction for interactor A = " + protein_a.uniprotAC + 
                                      " /interactor B = " + protein_b.uniprotAC )
        
        # Look at an interaction from B to A
        b_to_a_list = sql_session.query( PPINetworkInteraction ).filter( PPINetworkInteraction.interactorA == protein_b,
                                                       PPINetworkInteraction.interactorB == protein_a).all()
        
        # Check the result: if a single interaction is found, return it. If several are found, raise an exception
        if b_to_a_list != None and len( b_to_a_list ) > 0:
            if len( b_to_a_list ) == 1:
                return b_to_a_list[0]
            else:
                raise RainetException( "PPINetworkInteraction.__init__ : Found several protein interaction for interactor A = " + protein_b.uniprotAC + 
                                      " /interactor B = " + protein_a.uniprotAC + " /detection method = " + self.detectionMethod + " /pubmedID = " + self.pubmedID )
                
        return None
    
    # #
    # Add the object to SQLAlchemy session
    def add_to_session( self ):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self )
