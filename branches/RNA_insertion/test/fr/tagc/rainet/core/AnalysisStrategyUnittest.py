
import unittest
import os
import pandas as pd
import glob

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

    # Constants with default paramters        
    TOTAL_RNAS = 125 #200 total of all biotypes, 125 with used biotypes
    TOTAL_PROTS = 200
    TOTAL_PRIS = 24 #54 total of all biotypes, 27 with used biotypes, 24 when removing peptide redundancy
    TOTAL_PRIS_LINC_FILT = 1
        
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
        optionManager.set_option(OptionConstants.OPTION_OUTPUT_FOLDER, "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/" )
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, OptionConstants.DEFAULT_INTERACTION_SCORE)
        optionManager.set_option(OptionConstants.OPTION_RNA_BIOTYPES, OptionConstants.DEFAULT_RNA_BIOTYPES)
        optionManager.set_option(OptionConstants.OPTION_GENCODE, OptionConstants.DEFAULT_GENCODE)
        optionManager.set_option(OptionConstants.OPTION_EXPRESSION_VALUE_CUTOFF, OptionConstants.DEFAULT_EXPRESSION_VALUE_CUTOFF)

        
        # Set the level of verbosity
        Logger.get_instance().set_level(OptionManager.get_instance().get_option(OptionConstants.OPTION_VERBOSITY))

        # Setting up SQL manager
        SQLManager.get_instance().set_DBpath(DB_PATH)
        self.sql_session = SQLManager.get_instance().get_session()

        # setting up internal test folder paths
        self.expectedFolder = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_expected/Report"
        self.outputFolder = OptionManager.get_instance().get_option(OptionConstants.OPTION_OUTPUT_FOLDER ) + "/Report/"

        # create instance of strategy    
        self.strategy = AnalysisStrategy()
        
        # report only written for selected tests
        self.strategy.writeReportFile = 0
        

    # #
    # Test running AnalysisStrategy with default parameters
    def test_default_params(self):

        print "| test_default_params | "
        
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

        self.assertTrue(len(RNAs) == AnalysisStrategyUnittest.TOTAL_RNAS, "asserting if number of objects retrieved is correct") 
        self.assertTrue(len(Prots) == AnalysisStrategyUnittest.TOTAL_PROTS, "asserting if number of objects retrieved is correct") 
        self.assertTrue(len(PRIs) == AnalysisStrategyUnittest.TOTAL_PRIS, "asserting if number of objects retrieved is correct") 
  
    
    # #
    # Test filtering for mRNAs
    def test_RNA_filter_one(self):
  
        print "| test_RNA_filter_one | "
        
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_RNA_BIOTYPES, "protein_coding")
          
        self.strategy.execute()
  
        mRNAs = DataManager.get_instance().get_data(AnalysisStrategy.RNA_FILTER_KW)
          
        self.assertTrue(len(mRNAs) == 82, "asserting if number of objects retrieved is correct") 
   
        for mRNA in mRNAs:       
            self.assertTrue(isinstance(mRNA, MRNA), "check if the mRNA is instance of MRNA table/class")
  
    # #
    # Test filtering for specific subtypes of lncRNAs
    def test_RNA_filter_two(self):
  
        print "| test_RNA_filter_two | "
        
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_RNA_BIOTYPES, "lincRNA")
        self.strategy.execute()
  
        lincRNAs = DataManager.get_instance().get_data(AnalysisStrategy.RNA_FILTER_KW)
          
        self.assertTrue(len(lincRNAs) == 9, "asserting if number of objects retrieved is correct") 
    
        for lincRNA in lincRNAs:       
            self.assertTrue(isinstance(lincRNA, LncRNA), "check if the lncRNA is instance of LncRNA table/class")
  
        response = self.sql_session.query( LncRNA).filter( LncRNA.transcriptBiotype == "lincRNA" ).all()
          
        self.assertTrue(len(lincRNAs) == len(response), "asserting if number of objects in lncRNA default (off) option is same as querying directly lncRNA table") 
    
  
    # #
    # Test gencode filtering
    def test_RNA_filter_four(self):
        
        print "| test_RNA_filter_four | "
   
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_GENCODE, "1")
        self.strategy.execute()
   
        RNAs = DataManager.get_instance().get_data(AnalysisStrategy.RNA_FILTER_KW)
                
        self.assertTrue(len(RNAs) == 70, "asserting if number of object with gencode is correct") 
  
  
    # # 
    # Test simple PRI filter
    def test_PRI_filter_one(self):
   
        print "| test_PRI_filter_one | "

        optionManager = OptionManager.get_instance()        
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, "10.0")
        self.strategy.execute()
   
        PRIs = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)
        
        # Regarding peptide redundancy filter, there is two of such cases in test database, one is misc_RNA (is always filtered) and the other is lincRNA, so only one peptide (interaction) is removed
        
        self.assertTrue(len(PRIs) == 9, "asserting if number of interactions above certain interaction score is correct") 
  
  
    # #
    # Test PRI and RNA filter together
    def test_PRI_filter_two(self):
   
        print "| test_PRI_filter_two | "
   
        optionManager = OptionManager.get_instance()        
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, "28")
        optionManager.set_option(OptionConstants.OPTION_RNA_BIOTYPES, "lincRNA")
        optionManager.set_option(OptionConstants.OPTION_GENCODE, 1)
        self.strategy.execute()
   
        PRIs = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)
                   
        self.assertTrue(len(PRIs) == AnalysisStrategyUnittest.TOTAL_PRIS_LINC_FILT, "asserting if PRIs are affected by RNA-level filters") 


    # #
    # Test PRI expression filter
    def test_PRI_filter_three(self):
   
        print "| test_PRI_filter_three | "
   
        optionManager = OptionManager.get_instance()        
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, "28")
        optionManager.set_option(OptionConstants.OPTION_RNA_BIOTYPES, "lincRNA")
        optionManager.set_option(OptionConstants.OPTION_GENCODE, 1)
        optionManager.set_option(OptionConstants.OPTION_EXPRESSION_VALUE_CUTOFF, 2)

        self.strategy.execute()

        PRIs = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)
                                      
        self.assertTrue(len(PRIs) == AnalysisStrategyUnittest.TOTAL_PRIS_LINC_FILT, "asserting if PRIs are affected by RNA-level filters") 
   
        optionManager.set_option(OptionConstants.OPTION_EXPRESSION_VALUE_CUTOFF, 200)
        self.strategy.execute()

        PRIs = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)

        self.assertTrue(len(PRIs) == 0, "asserting if PRI is filtered with expression filter ") 


#     # #
#     # Test function to create report files with default parameters
#     # @unittest.skip("skipping")
#     def test_report_one(self):
# 
#         print "| test_report_one | "
# 
#         self.strategy.execute()
# 
#         # assert report files on filtering steps, if before and after filter have the same values
#         for report in [ AnalysisStrategy.REPORT_RNA_NUMBERS, AnalysisStrategy.REPORT_INTERACTION_NUMBERS]:
#             with open( self.outputFolder + report, "r") as out:
#                 header = out.readline()
#                 lineOne = out.readline().strip().split("\t")
#                 lineTwo = out.readline().strip().split("\t")
#                 print (lineOne,lineTwo)
#                 self.assertTrue(lineOne[1:] == lineTwo[1:], "assert that values before and after filter are the same with default parameters")
# 
#         # TODO: test for expression data


    # #
    # Test function to create report files with filters
    def test_report_two(self):

        print "| test_report_two | "

        optionManager = OptionManager.get_instance()        
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, "28")
        optionManager.set_option(OptionConstants.OPTION_RNA_BIOTYPES, "lincRNA")
        optionManager.set_option(OptionConstants.OPTION_GENCODE, 1)

        self.strategy.execute()
        
        # RNA numbers report
        table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_RNA_NUMBERS)
                
        self.assertTrue( table["Total_Genes"][0] == 198, "assert if number of Genes before filter is correct")
        self.assertTrue( table["Total_Genes"][1] == 8, "assert if number of Genes after filter is correct")        

        # Interaction numbers report
        table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_INTERACTION_NUMBERS)

        
        self.assertTrue( table["Total_interactions"][0] == 54,
                          "assert if number of total PRI before filter matches same value as other test")
        
        self.assertTrue( table["Total_interactions"][1] == AnalysisStrategyUnittest.TOTAL_PRIS_LINC_FILT,
                          "assert if number of total PRI after filter matches same value as other test")

        # Interaction scores report
        table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_INTERACTION_SCORES_BIOTYPE, header = None, sep = ",", skip_blank_lines = True)
                
        self.assertTrue( table[0][5] == "lincRNA", "assert if lincRNAs in file as they should")
        self.assertTrue( table[1][5] == 28.16, "assert if last value of lincRNA score is correct")

        # Interaction partners report
        table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_INTERACTIONS_SCORE_MATRIX, header = None, sep = "\t", skip_blank_lines = True)
        self.assertTrue( float(table[1][1]) == 28.16, "assert interaction value is correct")

        # Interaction partners report
        table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_INTERACTION_PARTNERS_BIOTYPE, header = None, sep = ",", skip_blank_lines = True)

        self.assertTrue( table[0][5] == "lincRNA", "assert if lincRNAs in file as they should")
        self.assertTrue( table[1][5] == 1, "assert if number of protein partners is correct")

        # RNA expression data presence report
        table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_RNA_EXPRESSION_DATA_PRESENCE, header = None, sep = "\t", skip_blank_lines = True)
        #confirmed using SQLLite
        self.assertTrue( table[1][1] == "6","assert if number of lincRNAs with expression data is correct")
        
        # RNA expression report
        table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_RNA_EXPRESSION, header = None, sep = "\t", skip_blank_lines = True)    
        #2  ENST00000426044  LncRNA            lincRNA            0.89, calculated using SQLLite and excel
        self.assertTrue( table[3][2] == "0.89","asserting if average expression value for given RNA is correct")





    # #
    # Test if report files are correctly formatted with manual inspection
    # Also added this so that I'm sure to remember to produce test for every report file I create
    def test_report_files(self):
        
        print "| test_report_files | "
        
        self.strategy.execute()

        # list of report files created with AnalysisStrategy       
        reportConstants = glob.glob( self.outputFolder + "/*")
        
        # assert if all files are as they should by manual inspection
        for report in reportConstants:
            with open(report, "r") as out:                
                with open(self.expectedFolder + "/" + os.path.basename(report), "r") as exp:
                    self.assertTrue(out.read() == exp.read(), "assert if report file is correct, by expected content comparison" )

        # Update helper:
        # 1) comment out "tearDown" function
        # 2) run the "test_default_params" function only
        # 3) manually check if the files are ok
        # 4) copy the test results files to expected folder: 
        # cp -r /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/ /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_expected
        # 5) uncomment "tearDown" and rerun all tests


    def test_sweaving(self):
  
        print "| test_sweaving | "
 
        # overwrite switch to write report file
        self.strategy.writeReportFile = 1
                    
        self.strategy.execute()


#     def test_extra(self):
#         
#         self.strategy.execute()

    # #
    # Runs after each test
    def tearDown(self):
                 
        # Wipe output folder
        cmd = "rm %s/*" % self.outputFolder
        os.system(cmd)
          
      


