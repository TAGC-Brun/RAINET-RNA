
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger

from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from sqlalchemy.orm.base import instance_dict


# #
# This class describes a Gene models obtained from querying Ensembl BioMart, to make correspondence between transcripts coming from the same Gene
#
class Gene( Base ):
    __tablename__ = 'Gene'

    # The Ensembl Gene ID
    geneID = Column( String, primary_key = True )
    # The base RNA
    transcriptList = relationship( 'RNA' , backref="Gene")


    # #
    # The Gene constructor
    # 
    # @param gene_ID : string - The provided Gene ID  
    def __init__( self, gene_ID):
        
        self.geneID = gene_ID

    # #
    # Add a RNA transcript to the list
    #
    # @param rna : string - A RNA instance
    def add_rna( self, rna ):

        if rna != None:
            self.transcriptList.append( rna )
         
