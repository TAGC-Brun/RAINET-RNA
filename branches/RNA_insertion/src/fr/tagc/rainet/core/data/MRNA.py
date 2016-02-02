
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey

from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager


class MRNA( RNA ):
    
    __tablename__ = 'MRNA'
    
    transcriptID = Column( String, ForeignKey('RNA.transcriptID'), primary_key=True)

    __mapper_args__ = { 'polymorphic_identity':'MRNA'}



    def __init__( self ):
        
        pass

