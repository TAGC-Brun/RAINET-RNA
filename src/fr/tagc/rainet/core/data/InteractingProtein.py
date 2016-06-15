
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey

from fr.tagc.rainet.core.util.sql.Base import Base

from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException


class InteractingProtein( Base ):
    
    __tablename__ = 'InteractingProtein'
    
    uniprotAC = Column( String, ForeignKey('Protein.uniprotAC'), primary_key=True)

    def __init__( self, protein_id ):
        
        #=======================================================================
        # Query Protein uniprotAC in respective DataManager object
        #=======================================================================
 
        from fr.tagc.rainet.core.util.data.DataManager import DataManager

        allProteins = DataManager.get_instance().get_data(DataConstants.PROT_ALL_KW)
                
        if protein_id in allProteins:
            self.uniprotAC = protein_id
        else:
            Logger.get_instance().warning( "\n InteractingProtein.init : Protein ID not found:\t" + str( protein_id) )
            raise NotRequiredInstantiationException( "InteractingProtein.init : no insertion needed, Protein object not found.")


    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
