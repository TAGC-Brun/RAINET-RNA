
from sqlalchemy import Column, String, Float, ForeignKey, ForeignKeyConstraint, PrimaryKeyConstraint
from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from fr.tagc.rainet.core.util.exception.RainetException import RainetException


# #
# This class describes a expression of transcript in a tissue
#
class RNATissueExpression( Base ):
    __tablename__ = 'RNATissueExpression'

    # The expressed RNA
    transcriptID = Column( String, ForeignKey( 'RNA.transcriptID', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)    
    # The tissue name
    tissueName = Column( String, ForeignKey( 'Tissue.tissueName', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)
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

        print (transcript_id)

        #=======================================================================        
        # Approach: read single file which contains the expression values and the tissues.
        # Create the tissue objects while reading the file and add the correspondence between
        # tissue and RNA in the Tissue object
        #=======================================================================
        
        sql_session = SQLManager.get_instance().get_session()
        
        #=======================================================================
        # Search for the RNA object
        #=======================================================================
        from fr.tagc.rainet.core.data.RNA import RNA

        rna_query = sql_session.query( RNA).filter( RNA.transcriptID == transcript_id).all()
         
        if rna_query != None and len( rna_query) > 0:
            if len( rna_query) == 1:
                rna = rna_query[0]
            else:
                raise RainetException( "RNATissueExpression.init : Abnormal number of RNAs found for transcriptID = " + transcript_id ) 
        else:
            raise RainetException( "RNATissueExpression.init : RNA not found for transcriptID = " + transcript_id ) 
         
        #=======================================================================
        # Build the Tissue objects related to the RNA expression
        # get instance of tissue and see if already present, if not, create new Gene entry
        #=======================================================================
        from fr.tagc.rainet.core.data.Tissue import Tissue

        tissue_query = sql_session.query( Tissue).filter( Tissue.tissueName == tissue_name).first()
                         
        # if no Tissue with that tissue name found, create one
        if tissue_query == None:
            tissue = Tissue ( tissue_name, source_db)
        else:
            tissue = tissue_query

        tissue.add_expressed_rna( rna)
        sql_session.add( tissue )
        
        #=======================================================================
        # Search for the inserted expression to add supplementary info
        #=======================================================================
        # --Check if a single RNATissueExpression is found. If not, raise an issue
        # --If a single RNATIssueExpression is found, add the complementary info

        expression_list = sql_session.query( RNATissueExpression).filter( RNATissueExpression.transcriptID == rna.transcriptID, RNATissueExpression.tissueName == tissue.tissueName).all()
         
        expression = None
        if expression_list != None and len( expression_list) > 0:
            if len( expression_list) == 1:
                expression = expression_list[0]
                try:
                    expression.expressionValue = float(expression_value)
                except ValueError as ve:
                    raise RainetException( "RNATissueExpression.__init__ : The expression value is not a float: " + str( expression_value ), ve )
                #sql_session.add( expression)
            else:
                raise RainetException( "RNATissueExpression.init : Abnormal number of RNATissueExpression found for rna= " + transcript_id + " and tissue=" + tissue_name + " : " + str( len( expression_list))) 
        else:
            raise RainetException( "RNATissueExpression.init : No RNATissueExpression found for rna = " + transcript_id + " and tissue=" + tissue_name)

        raise NotRequiredInstantiationException( "RNATissueExpression.init : RNATissueExpression objects do not have to be inserted by __init__ since they are created by Tissue to RNA association table.")
  
              
            
    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
    
