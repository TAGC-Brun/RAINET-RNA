
import unittest
import os
import numpy
import pandas as pd
import glob

from statsmodels.sandbox.stats.multicomp import multipletests

from fr.tagc.rainet.core.Rainet import Rainet
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.data.DataManager import DataManager
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.execution.AnalysisStrategy import AnalysisStrategy
from fr.tagc.rainet.core.execution.EnrichmentAnalysisStategy import EnrichmentAnalysisStrategy

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

#from UnittestConstants import *
from fr.tagc.rainet.core.util.exception.RainetException import RainetException

# #
# Unittesting the Rainet AnalysisStrategy on a specific validated database. 
#
# Before running these unittests one needs a Rainet database, populated with test data as in following command:
# Rainet.py Insertion -s human -d /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/rainet_testing_DB.sqlite -i /home/diogo/workspace/tagc-rainet-RNA/resources/insertion_human_rna_test.ini -f
#
class EnrichmentAnalysisStrategyUnittest(unittest.TestCase):

    # Constants with default paramters        
        
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
        optionManager.set_option(OptionConstants.OPTION_DB_NAME, "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/rainet_testing_DB.sqlite") #rainet2016-06-17.human_expression_wPRI.sqlite")
        optionManager.set_option(OptionConstants.OPTION_SPECIES, "human")
        optionManager.set_option(OptionConstants.OPTION_OUTPUT_FOLDER, "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/enrichmentAnalysis/" )
        optionManager.set_option(OptionConstants.OPTION_ANNOTATION_TABLE, "NetworkModule" )
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_PROTEIN_ANNOTATION, OptionConstants.DEFAULT_MINIMUM_PROTEIN_ANNOTATION)
        optionManager.set_option(OptionConstants.OPTION_MINIMUM_PROTEIN_INTERACTION, OptionConstants.DEFAULT_MINIMUM_PROTEIN_INTERACTION)
        optionManager.set_option(OptionConstants.OPTION_NUMBER_RANDOMIZATIONS, OptionConstants.DEFAULT_NUMBER_RANDOMIZATIONS)
        
        # Set the level of verbosity
        Logger.get_instance().set_level(OptionManager.get_instance().get_option(OptionConstants.OPTION_VERBOSITY))

        # Setting up SQL manager
        SQLManager.get_instance().set_DBpath(OptionManager.get_instance().get_option(OptionConstants.OPTION_DB_NAME))
        self.sql_session = SQLManager.get_instance().get_session()

        # setting up internal test folder paths
        self.expectedFolder = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/enrichmentAnalysis/test_expected/Report"
        self.outputFolder = OptionManager.get_instance().get_option(OptionConstants.OPTION_OUTPUT_FOLDER ) + "/Report/"

        # create instance of strategy    
        self.run = EnrichmentAnalysisStrategy()
                
                
    def test_default_params(self):
    
        print "| test_default_params | "
        
        self.run.execute()

  
    def test_multiple_correction(self):

        print "| test_multiple_correction | "

        # compare against R      
        pvalues = [0.01,0.01,0.01,0.01,0.4,0.0001] 
        correctedPvalues = list( self.run.multiple_test_correction(pvalues, meth = "fdr_bh"))

        self.assertTrue( correctedPvalues == [0.012, 0.012, 0.012, 0.012, 0.40000000000000002, 0.00060000000000000006], "corrected values are different from R p.adjust method")


    def test_count_sign_tests(self):
            
        print "| test_count_sign_tests | "
        
        nTests = 100
        
        tests = numpy.empty( nTests, object)
        
        # create fake tests
        # 3/4s have significant results, 1/2 are skipped
        hyperResults = [0] * (nTests / 4) + [0.01] * (nTests / 2) + [0.5] * (nTests / 4)
        skipTests = [0] * (nTests / 2) + [1] * (nTests / 2)
        
        for i in xrange(0, nTests):
            text = "%s\t%s\t%s\t%s\t%s\t%s\t%e" % ( "rna_id", "annotID", "x", "m", "k", skipTests[ i], hyperResults[ i] )
            tests[ i] = text.split("\t")

        testsCorrected = self.run.correct_pvalues( nTests, hyperResults, tests)

        self.assertTrue( len( testsCorrected) == len( tests), "number of tests should be the same regardless of correction")        
        self.assertTrue( len( testsCorrected[ 0]) -2 == len( tests[ 0]), "correct_pvalues should have added two extra columns to output"  )

        self.assertTrue( type( testsCorrected[ 0][ self.run.SIGN_COLUMN] == type( "qwerty")), "ensuring object type" )
        self.assertTrue( type( testsCorrected[ 0][ self.run.WARNING_COLUMN] == type( "qwerty")), "ensuring object type" )
        self.assertTrue( testsCorrected[ 0][ self.run.SIGN_COLUMN] == "1" or testsCorrected[ 0][ self.run.SIGN_COLUMN] == "0", "ensuring object type" )
        self.assertTrue( testsCorrected[ 0][ self.run.WARNING_COLUMN] == "1" or testsCorrected[ 0][ self.run.WARNING_COLUMN] == "0", "ensuring object type" )
        
        signResults = self.run.count_sign_tests( testsCorrected)      

        self.assertTrue( signResults[ 0] == 75, "assert correct number of significant results based on input")
        self.assertTrue( signResults[ 1] == 50, "assert correct number of significant results based on input")


    def test_randomize_annotation(self):
    
        annotDict = { "one" : ["P1","P2","P3"], "two" : ["P2","P4"], "three" : ["P5"]}
    
        randomAnnotDict = self.run.randomize_annotation( annotDict)

        print annotDict
        print randomAnnotDict



#     # #
#     # Runs after each test
#     def tearDown(self):
#                   
#         # Wipe output folder
#         cmd = "rm %s/*" % self.outputFolder
#         os.system(cmd)
          
      


