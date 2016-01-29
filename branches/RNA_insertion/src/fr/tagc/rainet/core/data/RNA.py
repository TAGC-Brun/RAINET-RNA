
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger

from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.data.Gene import Gene

from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from sqlalchemy.orm.base import instance_dict


# #
# This class describes a RNA models obtained from querying Ensembl BioMart
#
class RNA( Base ):
    __tablename__ = 'RNA'

    # The Transcript ID
    transcriptID = Column( String, primary_key = True )
    # The parent gene Gene ID
    geneID = Column ( String, ForeignKey("Gene.geneID") )
    # The Peptide ID
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
    transcriptGencodeBasic = Column ( Boolean )
    # The genomic start position coordinate of transcript
    transcriptStart = Column( Integer )
    # The genomic end position coordinate of transcript
    transcriptEnd = Column( Integer )
    # The genomic strandness of transcript
    transcriptStrand = Column ( Integer )
    # The chromosome name of transcript
    chromosomeName = Column ( String )
    # The percentage of GC content in transcript
    percentageGCContent = Column ( Float )

    # The list of RNA cross references
    crossReferences = relationship( 'RNACrossReference', backref = 'RNA' )
    # The list of Protein-RNA interactions
    pri = relationship("ProteinRNAInteraction")
    # The subtype of RNA (link to subtables)
    type = Column ( String )
 
    #For mapping to subtypes
    __mapper_args__ = {'polymorphic_identity':'RNA', 'polymorphic_on':type }


    # #
    # The RNA constructor.
    # Note: RNA class is not instantiated, instead, RNA subtypes (mRNA, lncRNA, otherRNA) are instantiated. 
    # Note: Gene table is filled here 
    #
    # @param transcript_ID: String - The Transcript ID (MANDATORY)
    # @param gene_ID: String - The parent gene Gene ID
    # @param peptide_ID: String - The Peptide ID
    # @param transcript_biotype: String - The type/category of the transcript (MANDATORY)
    # @param transcript_length: Integer - The length of the transcript
    # @param transcript_source: String - The annotation source of the transcript
    # @param transcript_status: String - The confidence status of the transcript
    # @param transcript_tsl: String - The Ensembl confidence level of the transcript
    # @param transcript_gencode_basic: Boolean - Presence or absence of this transcript on GENCODE Basic dataset
    # @param transcript_start: Integer - The genomic start position coordinate of transcript
    # @param transcript_end: Integer - The genomic end position coordinate of transcript
    # @param transcript_strand: Integer - The genomic strandness of transcript
    # @param chromosome_name: String - The chromosome name of transcript
    # @param percentage_GC_content: Float - The percentage of GC content in transcript
    def __init__( self, transcript_ID, gene_ID, peptide_ID, transcript_biotype, transcript_length, transcript_source, transcript_status, transcript_tsl, transcript_gencode_basic, transcript_start, transcript_end, transcript_strand, chromosome_name, percentage_GC_content):

        sql_session = SQLManager.get_instance().get_session()

        #=======================================================================
        # Choose and instantiate the appropriate RNA subtype
        #=======================================================================
        
        if transcript_biotype != "":
            if transcript_biotype in DataConstants.RNA_MRNA_BIOTYPE:
                from fr.tagc.rainet.core.data.MRNA import MRNA
                myRNA = MRNA()
            elif transcript_biotype in DataConstants.RNA_LNCRNA_BIOTYPE:
                from fr.tagc.rainet.core.data.LncRNA import LncRNA
                myRNA = LncRNA()
            else:
                from fr.tagc.rainet.core.data.OtherRNA import OtherRNA
                myRNA = OtherRNA()
        else:
            raise RainetException( "RNA.__init__ : The value of transcript biotype is empty: " + str( transcript_gencode_basic ) )   
        
        #=======================================================================
        # Fill the main RNA variables
        #=======================================================================
        
        myRNA.transcriptBiotype = transcript_biotype
        
        if transcript_ID != "":
            myRNA.transcriptID = transcript_ID
        else:
            raise RainetException( "RNA.__init__ : The value of transcript ID is empty: " + str( transcript_ID ) )   
        
        myRNA.peptideID = peptide_ID
        
        try:
            myRNA.transcriptLength = int( transcript_length )
        except ValueError as ve:
            raise RainetException( "RNA.__init__ : The value of transcript length is not an integer: " + str( transcript_length ), ve )
        
        myRNA.transcriptSource = transcript_source
        
        myRNA.transcriptStatus = transcript_status
        
        myRNA.transcriptTsl = transcript_tsl
        
        if transcript_gencode_basic == "GENCODE basic":
            myRNA.transcriptGencodeBasic = 1
        elif transcript_gencode_basic == "":
            myRNA.transcriptGencodeBasic = 0
        else:
            raise RainetException( "RNA.__init__ : The value of GENCODE basic is not identified: " + str( transcript_gencode_basic ))   
                 
        try:
            myRNA.transcriptStart = int( transcript_start )
        except ValueError as ve:
            raise RainetException( "RNA.__init__ : The value of transcript start is not an integer: " + str( transcript_start ), ve )
        
        try:
            myRNA.transcriptEnd = int( transcript_end )
        except ValueError as ve:
            raise RainetException( "RNA.__init__ : The value of transcript end is not an integer: " + str( transcript_end ), ve )
        
        try:
            myRNA.transcriptStrand = int( transcript_strand )
        except ValueError as ve:
            raise RainetException( "RNA.__init__ : The value of transcript strand is not an integer: " + str( transcript_strand ), ve )
        
        myRNA.chromosomeName = chromosome_name
        
        myRNA.percentageGCContent = percentage_GC_content
        
        try:
            myRNA.percentageGCContent = float ( percentage_GC_content )
        except ValueError as ve:
            raise RainetException( "RNA.__init__ : The value of GC content percentage end is not a float: " + str( percentage_GC_content ), ve )

        
        #=======================================================================
        # Build the Gene objects related to the RNA
        # each Gene can contain several transcripts. get instance of gene and see if already present, if not, create new Gene entry
        #=======================================================================

        gene = sql_session.query( Gene ).filter( Gene.geneID == gene_ID).first()

        if gene == None: #if no Gene with that Gene ID found, create one
            gene = Gene( gene_ID )
 
        gene.add_rna( myRNA )
        
        sql_session.add( gene )
        

        # Add the respective RNA subclass object, note that RNA itself is not instantiated
        sql_session.add( myRNA )
        raise NotRequiredInstantiationException( "RNA.__init__: RNA instance is not required." )

    # #
    # Add a RNACrossReference to the RNA cross reference list
    #
    # @param cross_reference : string - A RNACrossReference instance
    def add_cross_reference( self, cross_reference ):
        
        if cross_reference != None:
            self.crossReferences.append( cross_reference )


