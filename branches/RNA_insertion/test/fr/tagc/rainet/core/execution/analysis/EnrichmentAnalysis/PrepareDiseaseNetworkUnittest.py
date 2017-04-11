
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
        self.proteinShare = -1

        # folder containing expected output files
        self.expectedFolder = "test_expected/"
        
        self.run = PrepareDiseaseNetwork(self.rainetDB, self.commonDiseaseFile, self.outputFolder, self.proteinShare )
            

    # #
    def test_read_common_disease_file(self):

        print "| test_read_common_disease_file | "

        self.run.read_common_disease_file()

        ## Test networkInfo variable
        
        netInfo = self.run.networkInfo
            
        #cut -f2,7 manually_curated_common_disease.txt| tail -n +2 | sort -u | wc -l  66
        self.assertTrue( len( netInfo) == 66)
        
        #ENST00000522875    PVT1    P50539    'prostate;cancer'    prostate cancer    PROSTATE CANCER |  |    77|NetworkModule    Q96SQ5,P50151,Q96M91,Q5TAP6,Q8N6Y0,Q8IUQ4,P40222,P52292,Q6NYC8,Q9GZM8,A0A0C4DH19,P15622,Q70UQ0,I3L4J3,Q6XQN6,O15066,H3BTW2,Q9H741,O75386,P50539,Q9NQ69,Q96CN9,Q9HAT0,Q8NF50,Q93062,Q8TD31,Q01850,O76015,Q96GY3,Q86WW8,Q6PCT2,Q96MT8,Q08AG5,Q68D86,Q96BZ8,Q9H1X1,Q8N2I2,Q9P2K3,Q9NX04,P63000,Q9NS73,P08670,P35900,O75069,Q9BVG8,Q96GV9,Q86YD7,A8K8P3,Q92552,Q86YD1,Q6P1W5,Q9Y2J4,Q8ND90,Q9NYB9,Q6ZU52,Q8TCX5,Q8N0S2,O60763,P14136,Q9H257,Q13077,Q9Y3B7,Q9BZW7,P25786,Q9BYV2,O15015,P43356,Q96EV8,Q9H6J7,Q8IXL7,Q9NP87,O60341,O75604,P51116,P24863,Q7Z3B3,Q2KHM9,Q9H788,Q5JR59,Q96CS2,Q92620,Q8N8B7,Q8WWY3,P14373,O60293,Q3B820,Q08117,Q9UBB9,P07196,P13762,Q9P0W2,Q6A162,Q86Y26,P50053,Q8IYE0,Q8NA61,Q53EZ4,Q8N1B4,Q5T5P2,Q9UIA0,Q8N5R6,Q05516,O95671,P06753,Q9Y2D8,Q9H875,O94972,P19012,Q8IWZ5,P63218,Q14687,Q14D33,Q6NX49,Q9BWJ2,Q86Y33,P41219,Q5VU43,Q9BXY8,Q9Y6D9,O75558,Q9P126,Q8IYX8,Q04695,Q8IUH5,Q92917,Q7Z6G3,P55081,Q9P0T4,Q9ULM2,Q9BUI4,Q9UPU9,Q96PV4,O43482,P36406,P12036,O00560,Q08379,Q07002,Q13976,Q9UJV3,Q5GH77,Q8IYF3,Q8N3L3,Q9Y4E8,P98175,Q96EA4,Q15323,P25791,Q13137,Q15654,Q14134,Q14135,Q8N720,Q8WUT1,Q8IY82,P48668,Q99996,O95990,P02538,Q9H2F5,O43347,Q9BRK4,Q8NDB6,P08727,Q5W0B1,Q15834,Q9Y2I6,Q96HB5,Q2TBE0,Q96I25,Q6IPM2,Q8TD10,Q2TAC2,Q8WVF2,Q15561,P43360,Q9NU19,P43364,Q12933,Q9H0W8,Q14129,Q8NHQ1,Q16512,Q99750,Q92752,Q9Y247,Q9HBZ2    O75069,Q6NX49,Q5W0B1,Q8N6Y0,Q8N1B4,Q6IPM2,Q8N3L3,Q9Y6D9,Q9HBZ2,O15066,Q96MT8

        self.assertTrue( netInfo[ "PVT1_77|NetworkModule"][0] == "PVT1")
        self.assertTrue( netInfo[ "PVT1_77|NetworkModule"][1] == "77|NetworkModule")
        self.assertTrue( netInfo[ "PVT1_77|NetworkModule"][2] == str(11))
        self.assertTrue( netInfo[ "PVT1_77|NetworkModule"][4] == "prostate cancer")

        # PVT1 vs 350|NetworkModule has several disease associations

        self.assertTrue( len(netInfo[ "PVT1_350|NetworkModule"][4].split(",")) == 7)



    # #
    def test_calculate_protein_share(self):

        print "| test_calculate_protein_share | "

        self.run.proteinShare = 1

        proteinsPerComplex = { "comp1" : set( ["p1", "p2", "p3"]), "comp2" : set( ["p3", "p5"]), "comp3" : set( ["p6"])}

        ## Test complexShare variable
        
        complexShare = self.run.calculate_protein_share(proteinsPerComplex)
        
        self.assertTrue( len( complexShare) == 1)


        self.run.proteinShare = 0
        self.run.read_common_disease_file()
        self.assertTrue( len( self.run.complexShare) == 50 * 51 / 2 )
        self.assertTrue( self.run.complexShare["28|NetworkModule_773|NetworkModule"] == 0)
        self.assertTrue( self.run.complexShare["79|NetworkModule_104|NetworkModule"] == 26)

 
    def test_run(self):
 
        print "| test_run | "
         
        self.run.run()
         
 
#     # #
#     # Runs after each test
#     def tearDown(self):
#                                   
#         # Wipe output folder
#         cmd = "rm %s/*" % self.outputFolder
#         os.system(cmd)
#             
       


