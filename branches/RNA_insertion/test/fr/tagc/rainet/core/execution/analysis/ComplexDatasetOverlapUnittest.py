
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

from fr.tagc.rainet.core.execution.analysis.EnrichmentAnalysis.ComplexDatasetOverlap import ComplexDatasetOverlap

# #
# Unittesting. 
#
class ComplexDatasetOverlapUnittest(unittest.TestCase):

    # Constants with default paramters        
        
    # #
    # Runs before each test
    # name of this function needs forcely to be 'setUp'
    def setUp(self):

        # Set the options

        self.rainetDB = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_backup/RNA/rainet2016-12-28.human_lncRNA_cutoff50/rainet2016-12-28.human_lncRNA_cutoff50.sqlite"        
        self.outputFolder = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/analysis/test_output/complex_dataset_overlap/"
        self.useInteractingProteins = 1
        self.listDatasets = ComplexDatasetOverlap.DEFAULT_DATASET_LIST
        self.highOverlapStat = ComplexDatasetOverlap.DEFAULT_HIGH_OVERLAP
        
        self.run = ComplexDatasetOverlap( self.rainetDB, self.outputFolder, self.useInteractingProteins, self.listDatasets, self.highOverlapStat)

        
    def test_read_rainet_db(self):

        print "| test_read_rainet_db | "
        
        self.run.read_rainet_db( )
        
        datasetDict = self.run.datasetDict
        
        self.assertTrue( len( datasetDict) == 5)

        self.assertTrue( len( datasetDict["BioplexCluster"]) == 354 )
        
        self.assertTrue( "Q16539" in datasetDict["BioplexCluster"]["254"])
        

    def test_group_overlap(self):

        print "| test_group_overlap | "
        
        ## No overlap
        
        group1 = {"P1", "P2", "P3"}
        group2 = {"P4", "P5", "P6"}
        
        perc1, perc2 = ComplexDatasetOverlap.group_overlap( group1, group2)

        self.assertTrue( perc1 == 0.0)
        self.assertTrue( perc2 == 0.0)

        ## 1 overlap
        
        group1 = {"P1", "P2", "P3"}
        group2 = {"P3", "P4", "P5", "P6",}
        
        perc1, perc2 = ComplexDatasetOverlap.group_overlap( group1, group2)

        self.assertTrue( "%.2f" % perc1 == "33.33")
        self.assertTrue( "%.2f" % perc2 == "25.00")
        
        ## full overlap
        group1 = {"P1", "P2", "P3"}
        group2 = {"P1", "P2", "P3"}
        
        perc1, perc2 = ComplexDatasetOverlap.group_overlap( group1, group2)

        self.assertTrue( "%.2f" % perc1 == "100.00")


    def test_intra_dataset_overlap(self):

        print "| test_intra_dataset_overlap | "

        self.run.listDatasets = ["BioplexCluster"]

        self.run.read_rainet_db( )
        
        self.run.intra_dataset_overlap()
        
        self.assertTrue( len( self.run.intraDatasetOverlapResults[ "BioplexCluster"]) == len( self.run.datasetDict["BioplexCluster"]))


        # Test: O00303 is found in two complexes, but the other proteins of complex 102a are only found in that complex
        # "102a","O00303"
        # "102a","Q8IUZ0"
        # "102a","Q9H6J7"
        # "25b","O00303"        

        testingResults = self.run.intraDatasetOverlapResults[ "BioplexCluster"]["102a"] #e.g. BioplexCluster    102a    0.09    1    0

        spl = testingResults.strip().split( "\t")
                
        mean = spl[2]
        anyOverlap = spl[3]
        highOverlap = spl[4]
               
        # test mean overlap
        expectedVal = "%.2f" % (33.33 / len( self.run.datasetDict["BioplexCluster"]) ) 
        self.assertTrue( mean == expectedVal)
 
        # test any overlap
        self.assertTrue( anyOverlap == "1")
        
        # test high overlap
        self.assertTrue( str(highOverlap) == "0")
        

        ### testing the highOverlapStat parameter, decreasing it should find other overlaps as high
        self.run.highOverlapStat = 30

        self.run.intra_dataset_overlap()
        
        testingResults = self.run.intraDatasetOverlapResults[ "BioplexCluster"]["102a"] #e.g. BioplexCluster    102a    0.09    1    0
              
        highOverlap = testingResults.strip().split( "\t")[4]
        
        # test high overlap
        self.assertTrue( str(highOverlap) == "1")

        
    def test_inter_dataset_overlap(self):

        print "| test_inter_dataset_overlap | "

        self.run.listDatasets = ["BioplexCluster", "NetworkModule"]

        self.run.read_rainet_db( )
        
        self.run.inter_dataset_overlap()
        
        self.assertTrue( len( self.run.interDatasetOverlapResults[ "BioplexCluster|NetworkModule"]) == len( self.run.datasetDict["BioplexCluster"]) )       
        self.assertTrue( len( self.run.interDatasetOverlapResults[ "NetworkModule|BioplexCluster"]) == len( self.run.datasetDict["NetworkModule"]) )       

        # Testing example  # BioplexCluster|NetworkModule    97b     0.29    4       1
        # 
        # Bioplex search:
        # "97b","P28347"
        # "97b","Q14135"
        # 
        # Network module search:
        # "10","P28347"
        # "824","P28347"
        # 
        # "77","Q14135"
        # "444","Q14135"
        # "824","Q14135"
        # 
        # # there are 4 network modules with overlap, one of them 100%, others with 50% 
        
        testingResults = self.run.interDatasetOverlapResults[ "BioplexCluster|NetworkModule"]["97b"] 
        
        spl = testingResults.strip().split( "\t")
                
        mean = spl[2]
        anyOverlap = spl[3]
        highOverlap = spl[4]
               
        # test mean overlap
        expectedVal = "%.2f" % ( ( 50 + 50 + 50 + 100.0) / len( self.run.datasetDict["NetworkModule"]) ) 
        self.assertTrue( mean == expectedVal)
 
        # test any overlap
        self.assertTrue( anyOverlap == "4")
        
        # test high overlap
        self.assertTrue( str(highOverlap) == "1")
        
        
        
        
#     def test_run(self):
#  
#         print "| test_run | "        
#  
#         self.run.run()



#     def test_extra(self):
#   
#         self.run.catrapidFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/cutoff50/Corum_Havugimana/storedInteractions.tsv"
#         self.run.numberRandomizations = 1
#         self.run.topPartners = 50
#         self.run.run()


     
#     # #
#     # Runs after each test
#     def tearDown(self):
#                                 
#         # Wipe output folder
#         cmd = "rm %s/*" % self.outputFolder
#         os.system(cmd)
            
       


