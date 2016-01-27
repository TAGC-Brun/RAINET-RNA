
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey, PrimaryKeyConstraint
from sqlalchemy import or_
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger

from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from sqlalchemy.orm.base import instance_dict


# #
# This class describes an ID cross reference for RNA models, obtained from querying Ensembl BioMart
#
class RNACrossReference( Base ):
    __tablename__ = 'RNACrossReference'

    # The base RNA
    transcriptID = Column( Integer, ForeignKey( 'RNA.transcriptID' ) )
    # The database the information was obtained from
    sourceDB = Column( String )
    # The cross reference
    crossReferenceID = Column( String )
    # Define the Composite PrimaryKey
    __table_args__ = (
        PrimaryKeyConstraint('transcriptID', 'sourceDB', "crossReferenceID"),
    )


    # The constructor of the class
    #
    # @param RNA_acc : string - The ID of the RNA
    # @param db_source : string - The database name from which the RNA cross reference comes from
    # @param cross_reference : string - The cross reference ID
    def __init__(self, RNA_acc, db_source, cross_reference):
        
        # Get a SQLalchemy session
        sql_session = SQLManager.get_instance().get_session()
        
        # Retrieve the list of RNAs corresponding to the provided accession number
        from fr.tagc.rainet.core.data.RNA import RNA
        RNA_list = sql_session.query( RNA ).filter( or_( RNA.transcriptID == RNA_acc)).all() #this returns python object, RNA instance

        # If a single RNA exists with the given accession number, create the new RNACrossReference 
        # object with the right value and and it to the RNA cross reference list
        if RNA_list != None and len( RNA_list) == 1 :
            RNA = RNA_list[0]
            self.sourceDB = db_source
            self.crossReferenceID = cross_reference
            RNA.add_cross_reference( self )
            sql_session.add( RNA )
        # If several RNAs are found with the accession number raise an issue
        # If no RNAs are found, raise a NotRequiredInstantiationException to indicate
        # the new RNACrossReference object does not have to be inserted in DB.
        else:
            if len( RNA_list) > 1:
                raise RainetException( "RNACrossReference.init : Abnormal number of RNAs found for accession number '" + RNA_acc + "' : " + str( len( RNA_list)) + " RNAs found.")
            else:
                raise NotRequiredInstantiationException( "RNACrossReference.init : No Corresponding RNA found in Database." )

    ##
    # Add the object to SQLAlchemy session if it is linked to a RNA
    def add_to_session(self):
    
        if self.transcriptID != None and self.transcriptID != "":
            sql_session = SQLManager.get_instance().get_session()
            sql_session.add( self)
    
