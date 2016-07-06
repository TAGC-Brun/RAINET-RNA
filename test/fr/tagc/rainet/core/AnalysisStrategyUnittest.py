
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
from fr.tagc.rainet.core.data.InteractingProtein import InteractingProtein
from fr.tagc.rainet.core.data.InteractingRNA import InteractingRNA
from fr.tagc.rainet.core.data.Tissue import Tissue

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
    TOTAL_PRIS = 12
    TOTAL_PRIS_LINC_FILT = 2
        
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
#        optionManager.set_option(OptionConstants.OPTION_DB_NAME, "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/rainet2016-06-17.human_expression_wPRI.sqlite")
        optionManager.set_option(OptionConstants.OPTION_DB_NAME, DB_PATH)
        optionManager.set_option(OptionConstants.OPTION_SPECIES, "human")
        optionManager.set_option(OptionConstants.OPTION_OUTPUT_FOLDER, "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/" )
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, OptionConstants.DEFAULT_INTERACTION_SCORE)
        optionManager.set_option(OptionConstants.OPTION_RNA_BIOTYPES, OptionConstants.DEFAULT_RNA_BIOTYPES)
        optionManager.set_option(OptionConstants.OPTION_GENCODE, OptionConstants.DEFAULT_GENCODE)
        optionManager.set_option(OptionConstants.OPTION_EXPRESSION_VALUE_CUTOFF, OptionConstants.DEFAULT_EXPRESSION_VALUE_CUTOFF)
        optionManager.set_option(OptionConstants.OPTION_EXPRESSION_TISSUE_CUTOFF, OptionConstants.DEFAULT_EXPRESSION_TISSUE_CUTOFF)

        
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

        # important to create new SQLManager session if changing database
        SQLManager.get_instance().close_session()

        optionManager = OptionManager.get_instance()        
#        optionManager.set_option(OptionConstants.OPTION_DB_NAME, "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/rainet2016-06-17.human_expression_wPRI.sqlite")
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, "28")
        optionManager.set_option(OptionConstants.OPTION_RNA_BIOTYPES, "lincRNA")
        optionManager.set_option(OptionConstants.OPTION_GENCODE, 1)
        self.strategy.execute()
   
        PRIs = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)
        
        print len(PRIs)
              
        self.assertTrue(len(PRIs) == AnalysisStrategyUnittest.TOTAL_PRIS_LINC_FILT, "asserting if PRIs are affected by RNA-level filters") 


    # #
    # Test PRI expression filter
    def test_PRI_filter_three(self):
   
        print "| test_PRI_filter_three | "
   
        # Overwrite default values
        optionManager = OptionManager.get_instance()        
        optionManager.set_option(OptionConstants.OPTION_DB_NAME, "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/rainet2016-06-17.human_expression_wPRI.sqlite")
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, "100")
        optionManager.set_option(OptionConstants.OPTION_RNA_BIOTYPES, "lincRNA")
        optionManager.set_option(OptionConstants.OPTION_GENCODE, 1)
        optionManager.set_option(OptionConstants.OPTION_EXPRESSION_VALUE_CUTOFF, 1.0)
        optionManager.set_option(OptionConstants.OPTION_EXPRESSION_TISSUE_CUTOFF, 1)

        # important to create new SQLManager session if changing database
        SQLManager.get_instance().close_session()

        # Run strategy step by step
        self.strategy = AnalysisStrategy()
        self.strategy.execute( run = 0 )
 
        self.strategy.filter_RNA()
        self.strategy.filter_protein()
        self.strategy.filter_PRI()
        
        selectedInteractions = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)
        self.assertTrue( len( selectedInteractions) == 4298, "assert if number of initial interactions is correct") 

        # run main function we want to test
        self.strategy.dump_filter_PRI_expression()
        
        # select count(distinct(proteinID)) from MRNA  --> 57366
        self.assertTrue( len( self.strategy.mRNADict) == 57366, "assert if number of mRNA to protein correspondence is correct" )

        # ENST00000005180 --> Q9Y258, ENST00000394905 --> Q9Y258 #confirmed in Ensembl82 website        
        self.assertTrue( set( self.strategy.mRNADict[ "Q9Y258"]) == set(["ENST00000005180", "ENST00000394905"]), "assert if a specific mRNA-protein correspondence is correct" )
        
        # select count(distinct(transcriptID)) from RNATIssueExpression --> 184278
        self.assertTrue( len( self.strategy.expressionDict) == 184278, "assert number of transcripts with expression data is correct")

        #
        tissues = [ str( tiss[0]) for tiss in self.sql_session.query( Tissue.tissueName ).all() ]
        self.assertTrue( len( self.strategy.expressionDict[ "ENST00000394905"]) == len( tissues), "assert number of expression values match number of tissues ")

        # grep ENST00000394905 transcript_expression_metrics_no_outliers.tsv
        # ENST00000394905    Skin - Sun Exposed (Lower leg)    0.124    0.199    0.000    1.60395658469    0.754
        boo = 0
        for tissTuple in self.strategy.expressionDict[ "ENST00000394905"]:
            if tissTuple == (0.124, "Skin - Sun Exposed (Lower leg)" ):
                boo = 1
        self.assertTrue( boo, "assert specific transcript expression data is correct")

        ## Really test expression filter

        # check protein expression of a specific protein
        # "ENST00000294652","Q5VWK0"
        # "ENST00000370040","Q5VWK0"
        # "ENST00000444143","Q5VWK0"
        # "ENST00000495380","Q5VWK0"
        self.assertTrue( len( self.strategy.ProtMRNATissueExpressions["Q5VWK0"]) == 4)
        self.assertTrue( self.strategy.ProtMRNATissueExpressions["Q5VWK0"]["ENST00000370040"]["Liver"] == 0)
        self.assertTrue( self.strategy.ProtMRNATissueExpressions["Q5VWK0"]["ENST00000495380"]["Testis"] == 3.363)

        # ENST00000495380, an MRNA of Q5VWK0 protein, has expression value > 1.0 RPKM in testis, given by one of the mRNAs
        # This protein is not present in any other tissue
        #"ENST00000294652","Pituitary","0.012"
        #"ENST00000444143","Testis","0.488"
        #"ENST00000495380","Testis","3.363"

        # InteractingRNA we want to test as positive interaction: ENST00000413466
        # "ENST00000413466","Testis","1.006"

        # InteractingRNA we want to test as negative interaction: ENST00000423943
        # "ENST00000423943","Testis","0.512"

        proteinExpressionTissues = DataManager.get_instance().get_data(AnalysisStrategy.PROT_TISSUES_KW)        
        rnaExpressionTissues = DataManager.get_instance().get_data(AnalysisStrategy.RNA_TISSUES_KW)
        expressedInteractionsTissues = DataManager.get_instance().get_data(AnalysisStrategy.PRI_TISSUES_KW)
        
        self.assertTrue( "Q5VWK0" in proteinExpressionTissues["Testis"])  # key -> tissue, value -> set of protein IDs
        self.assertTrue( "Q5VWK0" not in proteinExpressionTissues["Pancreas"])  # key -> tissue, value -> set of protein IDs
        
        self.assertTrue( "ENST00000413466" in rnaExpressionTissues["Testis"]) # key -> tissue, value -> set of tx IDs        
        self.assertTrue( "ENST00000413466|Q5VWK0" in expressedInteractionsTissues, "assert if interaction passes cutoffs") # key -> transcriptID|proteinID (pair), value -> set of tissues
        self.assertTrue( len( expressedInteractionsTissues["ENST00000413466|Q5VWK0"]) == 1, "assert if number of tissues passing cutoff is correct") # key -> transcriptID|proteinID (pair), value -> set of tissues

        self.assertTrue( "ENST00000423943" not in rnaExpressionTissues["Testis"])
        self.assertTrue( "ENST00000423943|Q5VWK0" not in expressedInteractionsTissues, "assert that interaction is not present") # key -> transcriptID|proteinID (pair), value -> set of tissues
             
        newSelectedInteractions = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)

        self.assertTrue( len( newSelectedInteractions) <= len( selectedInteractions), "expression filtered interactions should be equal or less than initial ones")


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


#     # #
#     # Test function to create report files with filters
#     def test_report_two(self):
# 
#         print "| test_report_two | "
# 
#         optionManager = OptionManager.get_instance()        
#         optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, "28")
#         optionManager.set_option(OptionConstants.OPTION_RNA_BIOTYPES, "lincRNA")
#         optionManager.set_option(OptionConstants.OPTION_GENCODE, 1)
# 
#         self.strategy.execute()
#         
#         # RNA numbers report
#         table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_RNA_NUMBERS)
#                 
#         self.assertTrue( table["Total_Genes"][0] == 198, "assert if number of Genes before filter is correct")
#         self.assertTrue( table["Total_Genes"][1] == 8, "assert if number of Genes after filter is correct")        
# 
#         # Interaction numbers report
#         table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_INTERACTION_NUMBERS)
# 
#         
#         self.assertTrue( table["Total_interactions"][0] == 54,
#                           "assert if number of total PRI before filter matches same value as other test")
#         
#         self.assertTrue( table["Total_interactions"][1] == AnalysisStrategyUnittest.TOTAL_PRIS_LINC_FILT,
#                           "assert if number of total PRI after filter matches same value as other test")
# 
#         # Interaction scores report
#         table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_INTERACTION_SCORES_BIOTYPE, header = None, sep = ",", skip_blank_lines = True)
#                 
#         self.assertTrue( table[0][5] == "lincRNA", "assert if lincRNAs in file as they should")
#         self.assertTrue( table[1][5] == 28.16, "assert if last value of lincRNA score is correct")
# 
#         # Interaction partners report
#         table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_INTERACTIONS_SCORE_MATRIX, header = None, sep = "\t", skip_blank_lines = True)
#         self.assertTrue( float(table[1][1]) == 28.16, "assert interaction value is correct")
# 
#         # Interaction partners report
#         table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_INTERACTION_PARTNERS_BIOTYPE, header = None, sep = ",", skip_blank_lines = True)
# 
#         self.assertTrue( table[0][5] == "lincRNA", "assert if lincRNAs in file as they should")
#         self.assertTrue( table[1][5] == 1, "assert if number of protein partners is correct")
# 
#         # RNA expression data presence report
#         table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_RNA_EXPRESSION_DATA_PRESENCE, header = None, sep = "\t", skip_blank_lines = True)
#         #confirmed using SQLLite
#         self.assertTrue( table[1][1] == "6","assert if number of lincRNAs with expression data is correct")
#         
#         # RNA expression report
#         table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_RNA_EXPRESSION, header = None, sep = "\t", skip_blank_lines = True)    
#         #2  ENST00000426044  LncRNA            lincRNA            0.89, calculated using SQLLite and excel
#         self.assertTrue( table[3][2] == "0.89","asserting if average expression value for given RNA is correct")





#     # #
#     # Test if report files are correctly formatted with manual inspection
#     # Also added this so that I'm sure to remember to produce test for every report file I create
#     def test_report_files(self):
#         
#         print "| test_report_files | "
#         
#         self.strategy.execute()
# 
#         # list of report files created with AnalysisStrategy       
#         reportConstants = glob.glob( self.outputFolder + "/*")
#         
#         # assert if all files are as they should by manual inspection
#         for report in reportConstants:
#             with open(report, "r") as out:                
#                 with open(self.expectedFolder + "/" + os.path.basename(report), "r") as exp:
#                     self.assertTrue(out.read() == exp.read(), "assert if report file is correct, by expected content comparison" )
# 
#         # Update helper:
#         # 1) comment out "tearDown" function
#         # 2) run the "test_default_params" function only
#         # 3) manually check if the files are ok
#         # 4) copy the test results files to expected folder: 
#         # cp -r /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/ /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_expected
#         # 5) uncomment "tearDown" and rerun all tests
# 
# 
#     def test_sweaving(self):
#   
#         print "| test_sweaving | "
#  
#         # overwrite switch to write report file
#         self.strategy.writeReportFile = 1
#                     
#         self.strategy.execute()


#     # #
#     # Runs after each test
#     def tearDown(self):
#                   
#         # Wipe output folder
#         cmd = "rm %s/*" % self.outputFolder
#         os.system(cmd)
          
      


