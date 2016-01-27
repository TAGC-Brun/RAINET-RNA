
from sqlalchemy import Column, String, Integer, Boolean, Float
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger

from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from sqlalchemy.orm.base import instance_dict



# #
# This class describes a RNA models obtained from querying Ensembl BioMart
#
#
class RNA( Base ):
    __tablename__ = 'RNA'

    # The Ensembl Transcript ID
    rnaID = Column( String, primary_key = True )
    # The parent gene Ensembl Gene ID
#    parentGene = relationship( 'Gene', backref = 'RNA' ) #dr: this will be connected to a gene table, so will be foreign key
    # The Ensembl Peptide ID, if existing
    peptideID = Column( String )
    # The type/category of the transcript
    transcriptBiotype = Column( String )
    # The length of the transcript
    transcriptLength = Column( Integer )
    # The annotation source of the transcript
    transcriptSource = Column( String )
    # The confidence status of the transcript
    transcriptStatus = Column( String )
    # The Ensembl confidence level of the transcript
    transcriptTsl = Column ( String )
    # Presence or absence of this transcript on GENCODE Basic dataset
    transcriptGencodeBasic = Column ( String ) # dr: could change this to boolean
    # The genomic start position coordinate of transcript
    transcriptStart = Column( Integer )
    # The genomic end position coordinate of transcript
    transcriptEnd = Column( Integer )
    # The genomic strandness of transcript
    transcriptStrand = Column ( String ) # dr: could change this to boolean
    # The chromosome name of transcript
    chromosomeName = Column ( String )
    # The percentage of GC content in transcript
    percentageGCContent = Column ( Float )
 


    # #
    # The RNA constructor
    # 
    # @param 
    #### dr: NEED TO FILL documentation    
    def __init__( self, rna_ID, parent_gene, peptide_ID, transcript_biotype, transcript_length, transcript_source, transcript_status, transcript_tsl, transcript_gencode_basic, transcript_start, transcript_end, transcript_strand, chromosome_name, percentage_GC_content):
        
        #=======================================================================
        # Fill the main protein variables
        #=======================================================================

        self.rnaID = rna_ID
#        self.parentGene  = parent_gene
        self.peptideID = peptide_ID
        self.transcriptBiotype = transcript_biotype
        try:
            self.transcriptLength = int( transcript_length )
        except ValueError as ve:
            raise RainetException( "RNA.__init__ : The value of transcript length is not an integer: " + str( transcript_length ), ve )
        self.transcriptSource = transcript_source
        self.transcriptStatus = transcript_status
        self.transcriptTsl = transcript_tsl
        self.transcriptGencodeBasic = transcript_gencode_basic
        try:
            self.transcriptLength = int( transcript_start )
        except ValueError as ve:
            raise RainetException( "RNA.__init__ : The value of transcript start is not an integer: " + str( transcript_start ), ve )
        try:
            self.transcriptLength = int( transcript_end )
        except ValueError as ve:
            raise RainetException( "RNA.__init__ : The value of transcript end is not an integer: " + str( transcript_end ), ve )
        self.strand = transcript_strand
        self.chromosomeName = chromosome_name
        self.percentageGCContent = percentage_GC_content
        try:
            self.transcriptLength = float ( percentage_GC_content )
        except ValueError as ve:
            raise RainetException( "RNA.__init__ : The value of GC content percentage end is not a float: " + str( percentage_GC_content ), ve )

    # #
    # Add the object to SQLAlchemy session
    def add_to_session( self ):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self )

