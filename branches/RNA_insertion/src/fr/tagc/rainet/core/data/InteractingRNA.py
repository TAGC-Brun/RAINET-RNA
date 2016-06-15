
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey

from fr.tagc.rainet.core.util.sql.Base import Base

from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException


class InteractingRNA( Base ):
    
    __tablename__ = 'InteractingRNA'
    
    transcriptID = Column( String, ForeignKey('RNA.transcriptID'), primary_key=True)

    def __init__( self, transcript_id ):
        
        #=======================================================================
        # Query RNA transcript ID in respective DataManager object
        #=======================================================================
 
        from fr.tagc.rainet.core.util.data.DataManager import DataManager

        allRNAs = DataManager.get_instance().get_data(DataConstants.RNA_ALL_KW)
                
        if transcript_id in allRNAs:
            self.transcriptID = transcript_id
        else:
            Logger.get_instance().warning( "\n InteractingRNA.init : Transcript ID not found:\t" + str( transcript_id) )
            raise NotRequiredInstantiationException( "InteractingRNA.init : no insertion needed, RNA object not found.")


    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
