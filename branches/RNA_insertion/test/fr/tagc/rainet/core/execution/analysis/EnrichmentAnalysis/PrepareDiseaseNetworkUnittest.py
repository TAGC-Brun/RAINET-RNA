
import unittest
import os
import pandas as pd
import glob

from fr.tagc.rainet.core.execution.analysis.EnrichmentAnalysis.PrepareDiseaseNetwork import PrepareDiseaseNetwork

# #
# Unittesting 
#
class PrepareDiseaseNetworkUnittest(unittest.TestCase):
        
    # #
    # Runs before each test
    # name of this function needs forcely to be 'setUp'
    def setUp(self):

        # Set the options
        self.rainetDB = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_backup/RNA/rainet2016-12-28.human_lncRNA_cutoff50/rainet2016-12-28.human_lncRNA_cutoff50.sqlite"
        self.commonDiseaseFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_input/EnrichmentAnalysis/PrepareDiseaseNetwork/manually_curated_common_disease.txt"
        self.outputFolder = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/enrichmentAnalysis/PrepareDiseaseNetwork"

        # folder containing expected output files
        self.expectedFolder = "test_expected/"
        
        self.run = PrepareDiseaseNetwork(self.rainetDB, self.commonDiseaseFile, self.outputFolder )
            

    # #
    def test_read_common_disease_file(self):

        print "| test_read_common_disease_file | "

        self.run.read_common_disease_file()



# 
#     def test_run(self):
# 
#         print "| test_run | "
#         
#         self.run.run()
#         
# 
#     # #
#     # Runs after each test
#     def tearDown(self):
#                                   
#         # Wipe output folder
#         cmd = "rm %s/*" % self.outputFolder
#         os.system(cmd)
#             
       


