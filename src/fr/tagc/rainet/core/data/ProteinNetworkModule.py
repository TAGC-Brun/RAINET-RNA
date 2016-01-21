
from sqlalchemy import Column, String, ForeignKey

from fr.tagc.rainet.core.util.sql.Base import Base



# #
# This class describes the owning of a Protein to a NetworkModule
#
class ProteinNetworkModule( Base ):
    __tablename__ = 'ProteinNetworkModule'
    
    # The NetworkModule 
    networkModule_id = Column( String, ForeignKey( 'NetworkModule.moduleID', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)
    # The associated protein
    protein_id = Column( String, ForeignKey( 'Protein.uniprotAC', onupdate="CASCADE", ondelete="CASCADE"), primary_key=True)
