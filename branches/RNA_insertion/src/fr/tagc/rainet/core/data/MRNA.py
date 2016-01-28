
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey

from fr.tagc.rainet.core.data.RNA import RNA

class MRNA( RNA ):
    
    __tablename__ = 'MRNA'
    
    transcriptID = Column( String, ForeignKey('RNA.transcriptID'), primary_key=True)
    
    __mapper_args__ = { 'polymorphic_identity':'MRNA'}

    def __init__( self ):
        pass

