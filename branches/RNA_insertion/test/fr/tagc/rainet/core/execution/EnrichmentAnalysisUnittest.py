
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
        self.run.sql_session = SQLManager.get_instance().get_session()
                

                
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
            
        print "| test_randomize_annotation | "
        
        annotDict = { "one" : ["P1","P2","P3"], "two" : ["P2","P4"], "three" : ["P5"], "four" : ["P1","P5"]}

        listOfProteins = [ prot for annot in annotDict for prot in annotDict[ annot]]
        lengths1Set = [ len( set( annotDict[ annot])) for annot in sorted(annotDict)]
        lengths1 = [ len( annotDict[ annot]) for annot in sorted(annotDict)]

        for i in xrange(1000):
            randomAnnotDict = self.run.randomize_annotation( annotDict, listOfProteins)
     
            lengths2 = [ len( randomAnnotDict[ annot]) for annot in sorted(randomAnnotDict)]
     
            self.assertTrue( lengths1 == lengths2, "assert that topology of annotation dictionary is unchanged with the shuffle")
    
#             lengths2Set = [ len( set( randomAnnotDict[ annot])) for annot in sorted(randomAnnotDict)]
#               
#             self.assertTrue( lengths1 == lengths1Set, "assert that there are no duplicate IDs in shuffled list")
#             self.assertTrue( lengths2 == lengths2Set, "assert that there are no duplicate IDs in shuffled list")

            listOfProteins2 = [ prot for annot in randomAnnotDict for prot in randomAnnotDict[ annot]]

            print sorted(set(listOfProteins))
            print sorted(set(listOfProteins2))

            self.assertTrue( set( listOfProteins) == set(listOfProteins2), "assert that there are no duplicate IDs in shuffled list")


    def test_empirical_pvalue(self):

        print "| test_empirical_pvalue | "
        
        nTests = 100

        observedSignificant = 10
        
        # random are 4/4 above
        listRandomSignificant = [11] * nTests
        assert len( listRandomSignificant) == nTests

        res, count = self.run.empirical_pvalue(listRandomSignificant, observedSignificant)

        self.assertTrue( res == 1.0, "assert that empirical pvalue is correct according to input")
        self.assertTrue( count == 0, "assert that empirical count above is correct according to input")

        # random are 0/4 above
        listRandomSignificant = [0] * nTests
        assert len( listRandomSignificant) == nTests

        res, count = self.run.empirical_pvalue(listRandomSignificant, observedSignificant)

        self.assertTrue( res == 0.0, "assert that empirical pvalue is correct according to input")
        self.assertTrue( count == 100, "assert that empirical count above is correct according to input")

        # random are 1/4 below, 1/4 with same, 2/4 above
        listRandomSignificant = [0] * (nTests / 4) + [10] * (nTests / 4) + [50] * (nTests / 2)
        assert len( listRandomSignificant) == nTests

        res, count = self.run.empirical_pvalue(listRandomSignificant, observedSignificant)
        
        self.assertTrue( res == 0.75, "assert that empirical pvalue is correct according to input")
        self.assertTrue( count == 25, "assert that empirical count above is correct according to input")


    def test_hypergeometric_test(self):
        
        print "| test_hypergeometric_test | "

        x, m, n, k = 5, 10, 50, 10

        res = "%.5f" % self.run.hypergeometric_test( x, m, n, k)

        self.assertTrue( res == "0.00067")


    def test_get_interaction_data(self):
        
        print "| test_get_interaction_data | "

        self.run.get_interaction_data()

        # store into accessible variables
        interactingProteins = DataManager.get_instance().get_data(EnrichmentAnalysisStrategy.PRI_PROT_KW)
        interactingRNAs = DataManager.get_instance().get_data(EnrichmentAnalysisStrategy.PRI_RNA_KW)
        interactions = DataManager.get_instance().get_data(EnrichmentAnalysisStrategy.PRI_KW)

        self.assertTrue( len( interactions) == 30, "assert number of retrieved interactions match expected" )
        self.assertTrue( len( interactingProteins) == 48, "assert number of retrieved interactions match expected" )
        self.assertTrue( len( interactingRNAs) == 17, "assert number of retrieved interactions match expected" )


    def test_annotation_report(self):
        
        print "| test_annotation_report | "

        self.run.execute()

        # For NetworkModule

        # Background should be overlap between interactingProteins and Proteins with annotation
        #SELECT count(distinct( InteractingProtein.uniprotAC)) FROM InteractingProtein, ProteinNetworkModule WHERE  InteractingProtein.uniprotAC == ProteinNetworkModule.protein_id
        self.assertTrue( len( self.run.backgroundProteins) == 37, "assert number of background proteins is correctly calculated")
        
        #SELECT count(distinct(ProteinNetworkModule.protein_id)) FROM ProteinNetworkModule -> 42
        self.assertTrue( len( self.run.protAnnotDict) == 42, "assert number of proteins with annotations is correctly calculated")
                
        self.assertTrue( len( self.run.allProteinsWithInteractionData) == len( DataManager.get_instance().get_data(EnrichmentAnalysisStrategy.PRI_PROT_KW)), "assert number of proteins with interactions is correctly calculated")

        #SELECT count(distinct(ProteinNetworkModule.networkModule_id)) FROM ProteinNetworkModule -> 82. But this does not count interactions.
        pool = {prot for annot in self.run.annotWithInteractionDict for prot in self.run.annotWithInteractionDict[annot]}

        self.assertTrue( len( pool) == len( self.run.backgroundProteins), "confirm number of proteins with annotation and interaction")
        
        # For KEGG Pathway 
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_ANNOTATION_TABLE, "KEGGPathway" )
       
        self.run.execute()
        
        #SELECT count(distinct(ProteinKEGGAnnotation.protein_id)) FROM ProteinKEGGAnnotation -> 22
        self.assertTrue( len( self.run.protAnnotDict) == 22, "assert number of proteins with annotations is correctly calculated")
              
        self.assertTrue( len( self.run.allProteinsWithInteractionData) == len( DataManager.get_instance().get_data(EnrichmentAnalysisStrategy.PRI_PROT_KW)), "assert number of proteins with interactions is correctly calculated")
 
        pool = {prot for annot in self.run.annotWithInteractionDict for prot in self.run.annotWithInteractionDict[annot]}
 
        self.assertTrue( len( pool) == len( self.run.backgroundProteins), "confirm number of proteins with annotation and interaction")


#     # #
#     # Runs after each test
#     def tearDown(self):
#                   
#         # Wipe output folder
#         cmd = "rm %s/*" % self.outputFolder
#         os.system(cmd)
          
      


