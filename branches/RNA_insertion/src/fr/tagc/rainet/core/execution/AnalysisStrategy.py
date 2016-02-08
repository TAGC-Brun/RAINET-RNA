
from sqlalchemy import or_,and_,distinct

from fr.tagc.rainet.core.execution.ExecutionStrategy import ExecutionStrategy
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.util import Constants as OptionConstantsConstants
from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.data.DataManager import DataManager

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
from fr.tagc.rainet.core.data.ProteinRNAInteractionCatRAPID import ProteinRNAInteractionCatRAPID

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
        self.species = OptionManager.get_instance().get_option(OptionConstants.OPTION_SPECIES)
        self.minimumInteractionScore =  OptionManager.get_instance().get_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE)
        self.transcriptBiotype = OptionManager.get_instance().get_option(OptionConstants.OPTION_TRANSCRIPT_BIOTYPE)
        self.lncRNABiotypes = OptionManager.get_instance().get_option(OptionConstants.OPTION_LNCRNA_BIOTYPES)


        # Check if minimumInteractionScore is float
        try:
            float(self.minimumInteractionScore)
        except TypeError:
            raise RainetException( "AnalysisStrategy.execute: Provided minimum interaction score is not a float.")

        # Check if transcript biotype is correct
        if self.transcriptBiotype not in OptionConstantsConstants.TRANSCRIPT_BIOTYPES:
            raise RainetException( "AnalysisStrategy.execute: Provided transcript biotype is not allowed: "+str(self.transcriptBiotype))

        # Check if list of lncRNA subtypes is correct
        self.lncRNABiotypes = self.lncRNABiotypes.split(",")
        for subtype in self.lncRNABiotypes:
            if subtype not in OptionConstantsConstants.LNCRNA_BIOTYPES:
                raise RainetException( "AnalysisStrategy.execute: Provided lncRNA biotype is not allowed: "+str(subtype))              


        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.DBPath)
        self.sql_session = SQLManager.get_instance().get_session()

        
        ### TO BE CONSTANTS
        # Filtering datasets
        self.RNAKeyword = "selectedRNAs"
        self.ProtKeyword = "selectedProteins"
        self.PRIKeyword = "selectedInteractions"
                
        self.analysis()
        
        
    def analysis(self):

        Logger.get_instance().info( "AnalysisStrategy.analysis: Starting..." )
        
        print ("ARGUMENTS:",self.DBPath,self.species,self.minimumInteractionScore,self.transcriptBiotype,self.lncRNABiotypes)

        self.filterRNA()
        
        self.filterProtein()
        
        self.filterPRI()

        self.enrichementAnalysis()


    # #
    # Filter RNA models
    def filterRNA(self):
        
        # Filter transcripts based on biotype

        #if LncRNA, loop over LncRNA table, if MRNA, loop through MRNA table
        #if antisense, lincRNA, loop through RNA table

        query = "query( "+self.transcriptBiotype+" ).all()"
        DataManager.get_instance().perform_query(self.RNAKeyword, query) 
        
        ## NOT SURE HOW TO DO THIS, PERFORM QUERY/FILTER ALL AT ONCE OR STORE DATA IN VARIABLE FOR SUBSEQUENT FITLERS?
        
#         if self.transcriptBiotype == "LncRNA":
#             #read online that we can query presence of IDs using a python list such as self.lncRNABiotypes
#             query = eval(query(self.transcriptBiotype).filter(self.transcriptBiotype.transcriptBiotype.in_(self.lncRNABiotypes).all() ) )
            

        # Filter transcripts based on gencode basic presence
        
        # Select one transcript isoform per gene

        pass

    # #
    # Filter protein models    
    def filterProtein(self):
        pass

    # #
    # Filter protein-RNA interactions
    def filterPRI(self):
       
        # Filter interactions based on minimumInteractionScore
        
        query = "query( ProteinRNAInteractionCatRAPID ).filter(ProteinRNAInteractionCatRAPID.interactionScore >= "+str(self.minimumInteractionScore)+").all()"

        DataManager.get_instance().perform_query(self.PRIKeyword, query) 

        
    def enrichementAnalysis(self):

        pass
#         selectedInteractions = DataManager.get_instance().get_data(self.PRIKeyword)
#         # Get protein IDs and transcript IDs present in filtered interactions
#         for inter in selectedInteractions:
#             print (inter.transcriptID,inter.proteinID)
#             
#             prot = self.sql_session.query( Protein ).filter( Protein.uniprotAC == inter.proteinID ).all()
# 
#             rna = self.sql_session.query( RNA ).filter( RNA.transcriptID == inter.transcriptID ).all()


