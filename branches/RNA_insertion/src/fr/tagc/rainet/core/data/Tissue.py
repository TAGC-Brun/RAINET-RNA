
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger

from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from fr.tagc.rainet.core.data.RNATissueExpression import RNATissueExpression
from sqlalchemy.sql.schema import PrimaryKeyConstraint, ForeignKeyConstraint


# #
# This class describes a the tissues used for expression data (e.g. obtained from GTEx project data), 
# to make correspondence between transcripts expression per tissue
#
class Tissue( Base ):
    __tablename__ = 'Tissue'

    # The tissue name (e.g. SMTS from GTEx)
    tissueName = Column( String, primary_key = True )

    # The database / dataset where the data comes from    
    sourceDB = Column( String, primary_key = True  )
    
    # #
    # The Tissue constructor
    # 
    # @param tissue_name : string - The provided tissue name 
    # @param source_db : the dataset / database where the data comes from
    def __init__( self, tissue_name, source_db):
        
        self.tissueName = tissue_name

        self.sourceDB = source_db

        self.add_to_session()
    
    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)



