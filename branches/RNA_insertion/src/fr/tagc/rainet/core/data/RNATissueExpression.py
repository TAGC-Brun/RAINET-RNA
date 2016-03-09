
from sqlalchemy import Column, String, Float, ForeignKey, ForeignKeyConstraint, PrimaryKeyConstraint
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.data.DataManager import DataManager
from fr.tagc.rainet.core.data import DataConstants

from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger


# #
# This class describes a expression of transcript in a tissue
#
class RNATissueExpression( Base ):
    __tablename__ = 'RNATissueExpression'

    # The expressed RNA
    transcriptID = Column( String, ForeignKey( 'RNA.transcriptID'), primary_key=True)    
    # The tissue name
    tissueName = Column( String, ForeignKey( 'Tissue.tissueName'), primary_key=True)
    # The expression value
    expressionValue = Column( Float )

    #
    # The constructor of the class
    #
    # @param tissue_name : string - the tissue name
    # @param transcript_id : string - the RNA transcript ID
    # @param expression_value : float - the actual expression value
    # @param source_db : string - database/dataset where the information comes from. This is given in the DataConstants.
    def __init__(self, transcript_id, tissue_name, expression_value, source_db):

        #=======================================================================        
        # Approach: read single file which contains the expression values and the tissues.
        # Create the tissue objects while reading the file and add the correspondence between
        # tissue and RNA in the Tissue
        #
        # Avoiding SQL queries for performance purposes
        #=======================================================================
        
        dt_manager = DataManager.get_instance()
        
        #=======================================================================
        # Search for the RNA object
        #=======================================================================

        # use previously created dictionary containing all RNA objects identified by transcriptID
        if transcript_id in dt_manager.get_data( DataConstants.RNA_ALL_KW):
            # Add RNA transcriptID correspondence to this instance
            self.transcriptID = transcript_id #here
        else:
            # Some transcript IDs from input expression file (e.g. GTEx) are deprecated in the Ensembl version used for the RNA models
            Logger.get_instance().warning( "RNATissueExpression.init : RNA not found for transcriptID = " + transcript_id ) 
            raise NotRequiredInstantiationException( "RNATissueExpression.init : RNATissueExpression objects not inserted since corresponding RNA is not found.")
         
        #=======================================================================
        # Build the Tissue objects related to the RNA expression
        # use temporary DataManager object to see if tissue is already present, if not, create new entry
        #=======================================================================
        from fr.tagc.rainet.core.data.Tissue import Tissue

        # Create data structure external to this instance to accumulate already processed tissue names
        kw = "tempSet"  
        if kw not in dt_manager.data.keys():
            # initialise data manager as a set
            dt_manager.store_data( kw, set())
  
        if tissue_name not in dt_manager.data[ kw]:
            Tissue( tissue_name, source_db)
            # add tissue names to data manager
            dt_manager.data[ kw].add( tissue_name)

        # Add tissue name correspondance to this instance
        self.tissueName = tissue_name
         
        #=======================================================================
        # Add supplementary info about the expression
        #=======================================================================

        try: 
            expression = float( expression_value)
        except ValueError as ve:
            raise RainetException( "RNATissueExpression.__init__ : The expression value is not a float: " + str( expression_value ), ve )

        self.expressionValue = expression
         
            
    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
    
