
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey

from fr.tagc.rainet.core.data.RNA import RNA

class LncRNA( RNA ):
    
    __tablename__ = 'LncRNA'
    
    transcriptID = Column( String, ForeignKey('RNA.transcriptID'), primary_key=True)
    
    __mapper_args__ = { 'polymorphic_identity':'LncRNA'}

    def __init__( self ):
        pass

