
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

from fr.tagc.rainet.core.execution.analysis.lncRNA_vs_mRNA.ProteinTargetRatio import ProteinTargetRatio

class ProteinTargetRatioUnittest(unittest.TestCase):
        
    # #
    # Runs before each test
    # name of this function needs forcely to be 'setUp'
    def setUp(self):

        # Set the options

        self.interactionsFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/mrna_vs_lncrna/test/storedInteractions_test.tsv"
        self.transcriptTypesFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/mrna_vs_lncrna/test/transcript_biotype.csv"      
        self.outputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/mrna_vs_lncrna/test/protein_target_ratio_test.out"
        
        self.run = ProteinTargetRatio( self.interactionsFile, self.transcriptTypesFile, self.outputFile)
        

    # #
    def test_read_transcript_types(self):

        print "| test_read_network_file | "

        self.run.read_transcript_types()

        print self.run.transcriptTypesFile
        
        
#        self.assertTrue( len(listOfNames) ==  61695)

        

     
#     # #
#     # Runs after each test
#     def tearDown(self):
#                                 
#         # Wipe output folder
#         cmd = "rm %s/*" % self.outputFolder
#         os.system(cmd)
            
       


