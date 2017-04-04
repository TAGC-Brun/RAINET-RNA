
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
        self.writeInteractions = 1
        self.batchSize = 1000000
        # self.extraMetrics = 1
        self.writeNormalisedInteractions = 0
        self.writeInteractionMatrix = 0
        self.booleanInteraction = 0
        self.sampleInteractions = -1

        # folder containing expected output files
        self.expectedFolder = "test_expected/"
        
        self.run = ReadCatrapid(self.catRAPIDFile, self.outputFolder, self.interactionCutoff, self.interactionFilterFile, 
                                self.rnaFilterFile, self.proteinFilterFile, self.writeInteractions, self.batchSize,
                                self.writeNormalisedInteractions, self.writeInteractionMatrix, self.booleanInteraction, self.sampleInteractions )
            

    # #
    def test_default_params(self):

        print "| test_default_params | "

        wantedPairs = self.run.read_interaction_filter_file( )

        proteinInteractionsMean, proteinInteractionsCounter = self.run.read_catrapid_file( wantedPairs, set(), set())
                
        # diogo@diogo-tower:~/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_input$ grep Q7Z419 * | grep ENST00000544329
        # catRAPID_interactions_test.txt:sp|Q7Z419|R144B_HUMAN ENST00000544329    17.96    0.47    0.00
        # filter_file.txt:Q7Z419    ENST00000544329

        self.assertTrue( round(proteinInteractionsMean["Q7Z419"],1) == 18.0)

        ## Check if mean is well calculated

        # catRAPID_interactions_test.txt:sp|Q7Z5L9|I2BP2_HUMAN ENST00000544329    27.14    0.69    0.01
        # filter_file.txt:Q7Z5L9    ENST00000544329
        # catRAPID_interactions_test.txt:sp|Q7Z5L9|I2BP2_HUMAN ENST00000544591    25.72    0.66    0.00
        # filter_file.txt:Q7Z5L9    ENST00000544591

        self.assertTrue( proteinInteractionsCounter["Q7Z5L9"] == 2)
        self.assertTrue( round(proteinInteractionsMean["Q7Z5L9"],1) == 26.4)


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
        self.assertTrue( round(proteinInteractionsMean[ "Q7Z569"], 1) == 258.3) 

        # check if output files are correct
        with open(self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME, "r") as out:                
            with open(self.expectedFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME, "r") as exp:
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

        with open(self.outputFolder + ReadCatrapid.RNA_INTERACTIONS_FILENAME, "r") as out:                
            with open(self.expectedFolder + "/extraMetrics" + ReadCatrapid.RNA_INTERACTIONS_FILENAME, "r") as exp:
                self.assertTrue(out.read() == exp.read(), "assert if report file is correct, by expected content comparison" )


        # cp /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_output/* /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_expected/extraMetrics


    # #
    def test_write_normalised_interactions(self):

        print "| test_write_normalised_interactions | "

        self.run.writeNormalisedInteractions = 1

        # no file filter
        wantedPairs = set()
        self.run.read_catrapid_file( wantedPairs, set(), set())

        self.run.write_normalised_interactions()

        with open(self.outputFolder + ReadCatrapid.NORMALISED_STORED_INTERACTIONS_FILENAME, "r") as out:                
            with open(self.expectedFolder + "/normalised/" + ReadCatrapid.NORMALISED_STORED_INTERACTIONS_FILENAME, "r") as exp:
                self.assertTrue(out.read() == exp.read(), "assert if report file is correct, by expected content comparison" )

        # cp /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_output/* /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_expected/normalised


    # #
    def test_min_max_normalisation(self):

        print "| test_min_max_normalisation | "

        normVal = self.run._min_max_normalisation(0, -10, 10)
     
        self.assertTrue( normVal == 0.5)


    # #
    def test_write_matrix_output(self):

        print "| test_write_matrix_output | "

        self.run.read_catrapid_file( set(), set(), set())

        self.run.write_matrix_output()
        
        with open(self.outputFolder + ReadCatrapid.INTERACTIONS_SCORE_MATRIX, "r") as out:
            count = 0
            for line in out:
                if count == 5:
                    break
                spl = line.split( "\t")
                self.assertTrue( len( spl) == 52, "asserting number of columns is number of proteins plus row header")
                
                if count == 0:
                    self.assertTrue( "Q7Z419" in spl[1], "asserting the sorting of file")
                
                if count == 1:
                    self.assertTrue( spl[0] == "ENST00000542804", "asserting the sorting of file")
                    # grep ENST00000542804 storedInteractions.tsv | grep Q7Z419
                    # sp|Q7Z419|R144B_HUMAN ENST00000542804    20.56    0.54    0.00

                    self.assertTrue( float( spl[1]) == 20.56, "asserting that score is correct")        
                count+=1
     

    # #
    def test_write_matrix_output_two(self):

        print "| test_write_matrix_output_two | "

        ## testing writing of matrix with 0 and 1 instead of interaction score

        # adding interaction cut off parameter
        self.run.booleanInteraction = 1

        self.run.read_catrapid_file( set(["Q7Z419_ENST00000542804", "Q7Z419_ENST00000542821", "Q7Z429_ENST00000542804"]), set(), set())

        self.run.write_matrix_output()
        
        with open(self.outputFolder + ReadCatrapid.INTERACTIONS_SCORE_MATRIX, "r") as out:
            count = 0
            for line in out:
                if count == 3:
                    break
                spl = line.split( "\t")
                            
                self.assertTrue( len( spl) == 3, "asserting number of columns is number of proteins plus row header")
                 
                if count == 0:
                    self.assertTrue( "Q7Z429" in spl[2], "asserting the sorting of file")
                 
                if count == 2:
                    self.assertTrue( spl[0] == "ENST00000542821", "asserting the sorting of file")
                    self.assertTrue( float( spl[2].strip()) == 0, "asserting that score is correct")        
                count+=1


    # #
    def test_params_three(self):

        print "| test_params_three | "

        self.run.sampleInteractions = 1 

        proteinInteractionsMean, proteinInteractionsCounter = self.run.read_catrapid_file( set(), set(), set())
        
        # if no other filters were applied, we should have 5 interactions for each item.
        itemCount = self.run.itemCount

        self.assertTrue( len( itemCount) == 2051)

        # ensure that at items have at least one interaction
        self.assertTrue( 0 not in itemCount.values())

        interactionSample = self.run.interactionSample

        # ensure that number of interactions taken is close to number of wanted number of interactions per item
        self.assertTrue( abs(2051 * self.run.sampleInteractions - len( interactionSample)) < 200 * self.run.sampleInteractions )

        # optionally I could use a random seed..

     
    # #
    # Runs after each test
    def tearDown(self):
                                 
        # Wipe output folder
        cmd = "rm %s/*" % self.outputFolder
        os.system(cmd)
            
       


