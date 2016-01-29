
from sqlalchemy import or_,and_,distinct

from fr.tagc.rainet.core.execution.ExecutionStrategy import ExecutionStrategy
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.DBParameter import DBParameter
from fr.tagc.rainet.core.data.GeneOntology import GeneOntology
from fr.tagc.rainet.core.data.GeneSymbol import GeneSymbol
from fr.tagc.rainet.core.data.KEGGPathway import KEGGPathway
from fr.tagc.rainet.core.data.NetworkModuleAnnotation import NetworkModuleAnnotation
from fr.tagc.rainet.core.data.NetworkModule import NetworkModule
from fr.tagc.rainet.core.data.OCGPartitionAnalysis import OCGPartitionAnalysis
from fr.tagc.rainet.core.data.PartitionAnalysis import PartitionAnalysis
from fr.tagc.rainet.core.data.PPINetworkInteraction import PPINetworkInteraction
from fr.tagc.rainet.core.data.PPINetwork import PPINetwork
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.ProteinDomain import ProteinDomain
from fr.tagc.rainet.core.data.ProteinGOAnnotation import ProteinGOAnnotation
from fr.tagc.rainet.core.data.ProteinInteraction import ProteinInteraction
from fr.tagc.rainet.core.data.ProteinIsoform import ProteinIsoform
from fr.tagc.rainet.core.data.ProteinKEGGAnnotation import ProteinKEGGAnnotation
from fr.tagc.rainet.core.data.ProteinNetworkModule import ProteinNetworkModule
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinReactomeAnnotation import ProteinReactomeAnnotation
from fr.tagc.rainet.core.data.ReactomePathway import ReactomePathway
from fr.tagc.rainet.core.data.SynonymGeneSymbol import SynonymGeneSymbol
from fr.tagc.rainet.core.data.Gene import Gene
from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.MRNA import MRNA
from fr.tagc.rainet.core.data.LncRNA import LncRNA
from fr.tagc.rainet.core.data.OtherRNA import OtherRNA
from fr.tagc.rainet.core.data.RNACrossReference import RNACrossReference
from fr.tagc.rainet.core.data.TableStatus import TableStatus
from fr.tagc.rainet.core.data import DataConstants


# #
# This class define the Strategy enabling user to perform DB query interactively
class AnalysisStrategy(ExecutionStrategy):

    # #
    # The Strategy execution method
    def execute(self):

        # Getting input arguments        
        self.DBPath = OptionManager.get_instance().get_option(OptionConstants.OPTION_DB_NAME)
        species = OptionManager.get_instance().get_option(OptionConstants.OPTION_SPECIES)

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.DBPath)
        sql_session = SQLManager.get_instance().get_session()
        
        # Check if the species in the DB correspond to the species given by the user
        SQLManager.check_species(sql_session, species, True)

        Logger.get_instance().info("AnalysisStrategy.execute : ...")  # dr: to write here


        # dr TESTING
        
        # Getting all transcripts of a Gene
        gene = sql_session.query( Gene ).filter( Gene.geneID == "ENSG00000004777").first()        
        txList = gene.transcriptList      
#         for tx in txList:
#             print ("%s\t%s\t%s\t%s" % (gene.geneID, tx.transcriptID,tx.transcriptBiotype,tx.transcriptLength) )

        # Looping over all lncRNAs
        lncRNAs = sql_session.query( LncRNA ).count()
        print ("Total number of lncRNAs:\t%i" % lncRNAs )
        distinctLncRNABiotypes = sql_session.query( LncRNA.transcriptBiotype ).distinct().all()
        for bioType in distinctLncRNABiotypes:
            bioType = str(bioType[0])
            count = sql_session.query( LncRNA ).filter( LncRNA.transcriptBiotype == bioType ).count()        
            print bioType,count

            
        # Looping over all mRNAs
        mRNAs = sql_session.query( MRNA ).count()
        print ("Total number of mRNAs:\t%i" % mRNAs )

        # Getting protein UniprotAC from ENSEMBL ID (ENSP) using Protein CrossReference table
        query = sql_session.query( ProteinCrossReference ).filter( ProteinCrossReference.sourceDB == "Ensembl" and ProteinCrossReference.crossReferenceID == "ENSP..").all()
        
        #        protein_list = sql_session.query( Protein).filter( or_( Protein.uniprotAC == protein_acc, Protein.uniprotID == protein_acc)).all()

#         self.analysis()
# 
# 
#     def analysis(self):
#         
#         print ("HERE")