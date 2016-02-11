
import unittest

from fr.tagc.rainet.core.Rainet import Rainet
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.data.DataManager import DataManager
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.execution.AnalysisStrategy import AnalysisStrategy

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

from UnittestConstants import *
from fr.tagc.rainet.core.util.exception.RainetException import RainetException

# #
# Unittesting the Rainet AnalysisStrategy on a specific validated database. 
#
# Before running these unittests one needs a Rainet database, populated with test data as in following command:
# Rainet.py Insertion -s human -d /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/rainet_testing_DB.sqlite -i /home/diogo/workspace/tagc-rainet-RNA/resources/insertion_human_rna_test.ini -f
#
class AnalysisStrategyUnittest(unittest.TestCase):
    
    # #
    # Runs before each test
    # name of this function needs forcely to be 'setUp'
    def setUp(self):
                
        # Set the options
        # Note: if running from command line / main script, the optionManager gets the default values, 
        # but in unittest we must set up all arguments, whether optional or not.
        # In the actual unittests I may override the default options, for testing.
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_VERBOSITY, "debug")
        optionManager.set_option(OptionConstants.OPTION_DB_NAME, DB_PATH)
        optionManager.set_option(OptionConstants.OPTION_SPECIES, "human")
        optionManager.set_option(OptionConstants.OPTION_TRANSCRIPT_BIOTYPE, OptionConstants.DEFAULT_BIOTYPE)
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, OptionConstants.DEFAULT_INTERACTION_SCORE)
        optionManager.set_option(OptionConstants.OPTION_LNCRNA_BIOTYPES, OptionConstants.DEFAULT_LNCRNA_BIOTYPES)
        optionManager.set_option(OptionConstants.OPTION_GENCODE, OptionConstants.DEFAULT_GENCODE)
        
        # Set the level of verbosity
        Logger.get_instance().set_level(OptionManager.get_instance().get_option(OptionConstants.OPTION_VERBOSITY))

        # Setting up SQL manager
        SQLManager.get_instance().set_DBpath(DB_PATH)
        self.sql_session = SQLManager.get_instance().get_session()

        # create instance of strategy    
        self.strategy = AnalysisStrategy()
                
        
    # #
    # Test running AnalysisStrategy with default parameters
    def test_default_params(self):

        #Note: if wanting to catch the exception (and its message) use the following (because unittest bypasses Rainet.py)
#         try:
#             self.strategy.execute()
#             DataManager.get_instance().perform_query("kww","query( faketable ).all()")
#         except RainetException as re:
#             Logger.get_instance().error(re.to_string())
#             self.assertTrue(False,"asserting if exception occurred during execution.")

        self.strategy.execute()

        RNAs = DataManager.get_instance().get_data(AnalysisStrategy.RNA_FILTER_KW)
        Prots = DataManager.get_instance().get_data(AnalysisStrategy.PROT_FILTER_KW)
        PRIs = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)

        self.assertTrue(len(RNAs) == 200, "asserting if number of objects retrieved is correct") 
        self.assertTrue(len(Prots) == 200, "asserting if number of objects retrieved is correct") 
        self.assertTrue(len(PRIs) == 54, "asserting if number of objects retrieved is correct") 
  
  
    # #
    # Test filtering for mRNAs
    def test_RNA_filter_one(self):
  
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_TRANSCRIPT_BIOTYPE, "MRNA")
          
        self.strategy.execute()
  
        mRNAs = DataManager.get_instance().get_data(AnalysisStrategy.RNA_FILTER_KW)
          
        self.assertTrue(len(mRNAs) == 82, "asserting if number of objects retrieved is correct") 
   
        for mRNA in mRNAs:       
            self.assertTrue(isinstance(mRNA, MRNA), "check if the mRNA is instance of MRNA table/class")
  
    # #
    # Test filtering for specific subtypes of lncRNAs
    def test_RNA_filter_two(self):
  
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_TRANSCRIPT_BIOTYPE, "LncRNA")
        optionManager.set_option(OptionConstants.OPTION_LNCRNA_BIOTYPES, "antisense,lincRNA")
        self.strategy.execute()
  
        lncRNAs = DataManager.get_instance().get_data(AnalysisStrategy.RNA_FILTER_KW)
          
        self.assertTrue(len(lncRNAs) == 18, "asserting if number of objects retrieved is correct") 
    
        for lncRNA in lncRNAs:       
            self.assertTrue(isinstance(lncRNA, LncRNA), "check if the lncRNA is instance of LncRNA table/class")
  
    # #
    # Test default lncRNA subtype option
    def test_RNA_filter_three(self):
  
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_TRANSCRIPT_BIOTYPE, "LncRNA")
        self.strategy.execute()
  
        lncRNAs = DataManager.get_instance().get_data(AnalysisStrategy.RNA_FILTER_KW)
        response = self.sql_session.query(LncRNA).all()
          
        self.assertTrue(len(lncRNAs) == len(response), "asserting if number of objects in lncRNA default (off) option is same as querying directly lncRNA table") 
    
  
    # #
    # Test gencode filtering
    def test_RNA_filter_four(self):
   
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_TRANSCRIPT_BIOTYPE, "LncRNA")
        optionManager.set_option(OptionConstants.OPTION_GENCODE, "1")
        self.strategy.execute()
   
        lncRNAs = DataManager.get_instance().get_data(AnalysisStrategy.RNA_FILTER_KW)
           
        self.assertTrue(len(lncRNAs) == 19, "asserting if number of objects in being both lncRNA and gencode is correct") 
  
  
    # # 
    # Test simple PRI filter
    def test_PRI_filter_one(self):
   
        optionManager = OptionManager.get_instance()        
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, "10.0")
        self.strategy.execute()
   
        PRIs = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)
           
        self.assertTrue(len(PRIs) == 18, "asserting if number of interactions above certain interaction score is correct") 
  
  
    # #
    # Test PRI and RNA filter together
    def test_PRI_filter_two(self):
   
        optionManager = OptionManager.get_instance()        
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, "28")
        optionManager.set_option(OptionConstants.OPTION_TRANSCRIPT_BIOTYPE, "LncRNA")
        optionManager.set_option(OptionConstants.OPTION_LNCRNA_BIOTYPES, "lincRNA")
        optionManager.set_option(OptionConstants.OPTION_GENCODE, 1)
        self.strategy.execute()
   
        PRIs = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)
   
        self.assertTrue(len(PRIs) == 2, "asserting if PRIs are affected by RNA-level filters") 


#     def test_extra(self):
#         
#         self.strategy.execute()
#         
#         self.strategy.check_isoforms()


    # #
    # Runs after each test
    def tearDown(self):
        pass
    
    
    
    
