
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

from fr.tagc.rainet.core.execution.analysis.NetworkScoreAnalysis.NetworkScoreAnalysis import NetworkScoreAnalysis


# #
# Unittesting the ReadCatrapid script on a specific validated dataset. 
#
class NetworkScoreAnalysisUnittest(unittest.TestCase):

    # Constants with default paramters        
        
    # #
    # Runs before each test
    # name of this function needs forcely to be 'setUp'
    def setUp(self):

        # Set the options

        self.networkFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/input_data/PROTEIN/human.binary.nr0.95.connected.noself.gr"
        self.catrapidFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/snrna/sn_expression_1.58_cutoff_15/storedInteractions.tsv"
        #self.catrapidFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/telomerase_plus_random_tx/telomerase_plus_random_tx_interactions.out"      
        self.rainetDBFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_backup/RNA/rainet2016-09-27.human_noPRI.sqlite"
        self.topPartners = 10
        self.outputFolder = "test_output/"
        self.numberRandomizations = 10

#         # folder containing expected output files
#         self.expectedFolder = "test_expected/"
        
        self.run = NetworkScoreAnalysis( self.networkFile,  self.catrapidFile,  self.rainetDBFile,  self.topPartners,  self.outputFolder,  self.numberRandomizations)
        

    # #
    def test_read_network_file(self):

        print "| test_read_network_file | "

        graph, listOfNames, listOfTuples, dictNames = self.run.read_network_file()
        
        # wc -l /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/input_data/PROTEIN/human.binary.nr0.95.connected.noself.gr
        # 61695
        
        self.assertTrue( len(listOfNames) ==  61695)

        # used comparelists pp1 pp2 all_combinations 0, checked lines of union.txt
        self.assertTrue( len(dictNames) ==  12318)
        

#     def test_protein_cross_references(self):
# 
#         print "| test_protein_cross_references | "
#         
#         self.run.protein_cross_references()
#         
#         proteinIDMappingDict = self.run.proteinIDMappingDict
# 
#         self.assertTrue( proteinIDMappingDict["1A23_HUMAN"] == "P30447")


    # #
    def test_calculate_protein_degree(self):

        print "| test_calculate_protein_degree | "

#        self.run.protein_cross_references()

        graph, listOfNames, listOfTuples, dictNames = self.run.read_network_file()
                     
        self.run.calculate_protein_degree()
   
        degreeDict = self.run.degreeDict

        # TCL1B_HUMAN = O95988                
        #         grep TCL1B  /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/input_data/PROTEIN/human.binary.nr0.95.connected.noself.gr 
        #         TCL1B_HUMAN    TNPO2_HUMAN
        #         AKT2_HUMAN    TCL1B_HUMAN

        self.assertTrue( degreeDict[ "TCL1B_HUMAN"] == 2)



    def test_read_catrapid_file(self):

        print "| test_read_catrapid_file | "        
                
        self.run.read_catrapid_file()

        rnaTargets = self.run.rnaTargets
        allProtSet = self.run.allProtSet       
        allRNASet = self.run.allRNASet

        # cut -f1 storedInteractions.tsv | cut -f1 -d" " | sort -u | wc -l
        # 4905

        self.assertTrue( len(allProtSet) == 4905)

        # cut -f1 storedInteractions.tsv | cut -f2 -d" " | sort -u | wc -l
        # 10

        self.assertTrue( len(allRNASet) == 10)

        self.assertTrue( "ENST00000383925" in rnaTargets)

        # grep ENST00000383925 storedInteractions.tsv | grep Q5T5D7
        # sp|Q5T5D7|ZN684_HUMAN ENST00000383925    17.75    0.47    0.00

        self.assertTrue( "ZN684_HUMAN" in rnaTargets["ENST00000383925"][17.75] )
        self.assertTrue( "E2AK2_HUMAN" in rnaTargets["ENST00000383925"][20.16] )


    def test_pick_top_proteins(self):

        print "| test_pick_top_proteins | "        

        self.run.read_network_file()

        self.run.read_catrapid_file()

        self.run.pick_top_proteins()

        rnaTops = self.run.rnaTops

        # Top 10 interactions for ENST00000384010, sorted in excel
        # sp|Q99661|KIF2C_HUMAN    ENST00000384010    51.05    0.95    0.38
        # sp|Q9H116|GZF1_HUMAN    ENST00000384010    50.14    0.95    0.35
        # sp|Q5TD94|RSH4A_HUMAN    ENST00000384010    49.52    0.94    0.29
        # sp|Q9Y6D9|MD1L1_HUMAN    ENST00000384010    45.69    0.92    0.18
        # sp|A0AV02|S12A8_HUMAN    ENST00000384010    45.44    0.92    0.18
        # sp|Q4G1C9|GRPL2_HUMAN    ENST00000384010    44.32    0.91    0.17
        # sp|Q53EV4|LRC23_HUMAN    ENST00000384010    43.98    0.91    0.14
        # sp|Q96AQ6|PBIP1_HUMAN    ENST00000384010    43.9    0.91    0.14
        # sp|Q70EK9|UBP51_HUMAN    ENST00000384010    43.68    0.91    0.14
        # sp|Q92541|RTF1_HUMAN    ENST00000384010    42.07    0.9    0.12

        self.assertTrue( rnaTops[ "ENST00000384010"] == ['KIF2C_HUMAN', 'MD1L1_HUMAN', 'GRPL2_HUMAN', 'PBIP1_HUMAN', 'RTF1_HUMAN', 'IGBP1_HUMAN', 'SHOT1_HUMAN', 'MUM1_HUMAN', 'CASC1_HUMAN', 'SSRP1_HUMAN'])


    def test_calculate_metric_for_rna(self):

        print "| test_calculate_metric_for_rna | "        

        ## changing parameters in other to be stable
        self.run.topPartners = 3
        distanceScore = {
                      1 : 100,
                      2 : 1.0,
                      3 : 0.01 }        
        NetworkScoreAnalysis.DISTANCE_SCORE = distanceScore

        # run
        self.run.read_network_file()
        self.run.calculate_protein_degree()
        self.run.read_catrapid_file()        
        self.run.pick_top_proteins()

        topProteins = self.run.rnaTops[ "ENST00000384382"]

        meanRNAShortestPath, lionelMetric = self.run._calculate_metric_for_rna( topProteins)

#         # top 3 proteins of ENST00000384382
#         MCM7_HUMAN
#         RPB4_HUMAN
#         PPM1G_HUMAN
        # In cytoscape I checked their first and second neighbours, then got network with nodes connecting the 3 initial proteins only.
        # Visually I can see that RPB4_HUMAN and MCM7_HUMAN are connected by one first neighbour (NEMO_HUMAN). 
        # NEMO_HUMAN is connected to RBM8A_HUMAN and MCM10_HUMAN as second neighbours
        # PPM1G_HUMAN is isolated in the subnetwork of second neighbours
        # In summary, for the sum of all proteins, the distances counts are: 1 -> 0 times, 2 -> 2 times, 3 -> 4 times.
                
        expectedResult = str(2.04)
        
        self.assertTrue( str( lionelMetric) == expectedResult )
        

    def test_get_sample_protein_degree(self):

        print "| test_get_sample_protein_degree | "               

        self.run.read_network_file()
        self.run.calculate_protein_degree()

        # picked some random proteins in the PPI
        topProteins = [ "GLP1R_HUMAN", "B2L12_HUMAN", "KAD5_HUMAN"]

        newProteinSet = self.run._get_sample_protein_degree( topProteins)

        originalDegrees = { self.run.degreeDict[ prot] for prot in topProteins}
        newDegrees = { self.run.degreeDict[ prot] for prot in newProteinSet}
        
        self.assertTrue( originalDegrees == newDegrees)

        newProteinSet2 = self.run._get_sample_protein_degree( topProteins)

        self.assertTrue( newProteinSet != newProteinSet2)


    def test_calculate_metrics(self):

        print "| test_calculate_metrics | "               

        self.run.read_network_file()
        self.run.calculate_protein_degree()
        self.run.read_catrapid_file()        
        self.run.pick_top_proteins()
        self.run.calculate_metrics()


    def test_run(self):

        print "| test_run | "        

        self.run.run()


     
    # #
    # Runs after each test
    def tearDown(self):
                                
        # Wipe output folder
        cmd = "rm %s/*" % self.outputFolder
        os.system(cmd)
            
       


