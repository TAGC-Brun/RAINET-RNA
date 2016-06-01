
import unittest
import os
import pandas as pd
import glob

from fr.tagc.rainet.core.execution.processing.catrapid.ReadCatrapid import ReadCatrapid

# #
# Unittesting the ReadCatrapid script on a specific validated dataset. 
#
class ReadCatrapidUnittest(unittest.TestCase):

    # Constants with default paramters        
        
    # #
    # Runs before each test
    # name of this function needs forcely to be 'setUp'
    def setUp(self):

        # Set the options

        self.catRAPIDFile = "test_input/catRAPID_interactions_test.txt"
        self.outputFolder = "test_output/"
        self.interactionCutoff = "OFF"
        self.interactionFilterFile = "test_input/filter_file.txt"
        self.rnaFilterFile = "test_input/rna_filter_file.txt"
        self.proteinFilterFile = "test_input/protein_filter_file.txt"
        self.batchSize = 1000000
        self.extraMetrics = 0

        # folder containing expected output files
        self.expectedFolder = "test_expected/"
        
        self.run = ReadCatrapid(self.catRAPIDFile, self.outputFolder, self.interactionCutoff, self.interactionFilterFile, self.rnaFilterFile, self.proteinFilterFile, self.batchSize, self.extraMetrics)
            

    # #
    def test_default_params(self):

        print "| test_default_params | "

        wantedPairs = self.run.read_interaction_filter_file( )

        proteinInteractionsMean, proteinInteractionsCounter = self.run.read_catrapid_file( wantedPairs, set(), set())
        
        # diogo@diogo-tower:~/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_input$ grep Q7Z419 * | grep ENST00000544329
        # catRAPID_interactions_test.txt:sp|Q7Z419|R144B_HUMAN ENST00000544329    17.96    0.47    0.00
        # filter_file.txt:Q7Z419    ENST00000544329

        self.assertTrue( proteinInteractionsMean["Q7Z419"] == 17.96)

        ## Check if mean is well calculated

        # catRAPID_interactions_test.txt:sp|Q7Z5L9|I2BP2_HUMAN ENST00000544329    27.14    0.69    0.01
        # filter_file.txt:Q7Z5L9    ENST00000544329
        # catRAPID_interactions_test.txt:sp|Q7Z5L9|I2BP2_HUMAN ENST00000544591    25.72    0.66    0.00
        # filter_file.txt:Q7Z5L9    ENST00000544591

        self.assertTrue( proteinInteractionsCounter["Q7Z5L9"] == 2)
        self.assertTrue( proteinInteractionsMean["Q7Z5L9"] == 26.43)


    # #
    def test_params_one(self):

        print "| test_params_one | "

        # no file filter
        wantedPairs = set()

        # no interaction cutoff        
        self.run.interactionCutoff = "OFF"

        proteinInteractionsMean, proteinInteractionsCounter = self.run.read_catrapid_file( wantedPairs, set(), set())

        # cut -f1 catRAPID_interactions_test.txt | cut -f2 -d"|" | sort -u | wc -l
        # 51

        self.assertTrue( len( proteinInteractionsMean) == 51)

        for prot in proteinInteractionsCounter:
            self.assertTrue( proteinInteractionsCounter[ prot] == 2000 or proteinInteractionsCounter[ prot] == 100)

    # #
    def test_params_two(self):

        print "| test_params_two | "

        # no file filter
        wantedPairs = set()

        # adding interaction cut off parameter
        self.run.interactionCutoff = 200

        proteinInteractionsMean, proteinInteractionsCounter = self.run.read_catrapid_file( wantedPairs, set(), set())

        self.assertTrue( "Q7Z419" not in proteinInteractionsCounter) 

        # grep Q7Z569 catRAPID_interactions_test.txt | cut -f2 | sort -un | l
        # there is 1 interactions above 200

        self.assertTrue( proteinInteractionsCounter[ "Q7Z569"] == 1) 
        self.assertTrue( proteinInteractionsMean[ "Q7Z569"] == 258.29) 

        # check if output files are correct
        with open(self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME + "1.tsv", "r") as out:                
            with open(self.expectedFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME + "1.tsv", "r") as exp:
                self.assertTrue(out.read() == exp.read(), "assert if report file is correct, by expected content comparison" )

        with open(self.outputFolder + ReadCatrapid.PROTEIN_INTERACTIONS_FILENAME, "r") as out:                
            with open(self.expectedFolder + ReadCatrapid.PROTEIN_INTERACTIONS_FILENAME, "r") as exp:
                self.assertTrue(out.read() == exp.read(), "assert if report file is correct, by expected content comparison" )

        # cp /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_output/* /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_expected/

    # #
    def test_read_interaction_filter_file(self):

        print "| test_read_interaction_filter_file | "

        wantedPairs = self.run.read_interaction_filter_file( )
        
        self.assertTrue( wantedPairs == set(['Q7Z419_ENST00000544329', 'Q7Z5L9_ENST00000544591', 'Q7Z5L9_ENST00000544329']), "asserting if object is correct") 
  

    # #
    def test_apply_filter_file(self):

        print "| test_apply_filter_file | "

        # RNA filter
        wantedRNAs = self.run.read_rna_filter_file( )
        
        proteinInteractionsMean, proteinInteractionsCounter = self.run.read_catrapid_file( set(), wantedRNAs, set())

        # grep ENST00000544089 catRAPID_interactions_test.txt | wc -l
        # 51
        # grep ENST00000547795 catRAPID_interactions_test.txt | wc -l
        # 50

        self.assertTrue( sum( proteinInteractionsCounter.values()) == 101, "asserting if correct number of total interactions")

        # Protein filter
        wantedProteins = self.run.read_protein_filter_file()
        
        proteinInteractionsMean, proteinInteractionsCounter = self.run.read_catrapid_file( set(), set(), wantedProteins)

        # grep Q7Z5L7 catRAPID_interactions_test.txt | wc -l
        # 2000
        # grep Q7Z429 catRAPID_interactions_test.txt | wc -l
        # 2000

        self.assertTrue( sum( proteinInteractionsCounter.values()) == 4000, "asserting if correct number of total interactions")

        # RNA and protein filter

        proteinInteractionsMean, proteinInteractionsCounter = self.run.read_catrapid_file( set(), wantedRNAs, wantedProteins)

        self.assertTrue( sum( proteinInteractionsCounter.values()) == 4, "asserting if correct number of total interactions")

        

    # #
    def test_extra_metrics(self):

        print "| test_extra_metrics | "

        # change tag to calculate extra metrics
        self.run.extraMetrics = 1

        proteinInteractionsMean, proteinInteractionsCounter = self.run.read_catrapid_file( set(), set(), set())

        with open(self.outputFolder + ReadCatrapid.PROTEIN_INTERACTIONS_FILENAME, "r") as out:                
            with open(self.expectedFolder + "/extraMetrics" + ReadCatrapid.PROTEIN_INTERACTIONS_FILENAME, "r") as exp:
                self.assertTrue(out.read() == exp.read(), "assert if report file is correct, by expected content comparison" )

        # cp /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_output/* /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_expected/extraMetrics

    
    # #
    # Runs after each test
    def tearDown(self):
                   
        # Wipe output folder
        cmd = "rm %s/*" % self.outputFolder
        os.system(cmd)
          
      


