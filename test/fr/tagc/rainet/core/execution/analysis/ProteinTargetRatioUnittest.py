
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

        print "| test_read_transcript_types | "

        self.run.read_transcript_types()

        self.assertTrue( len( self.run.transcriptType) == 216133)


    # #
    def test_read_interaction_file(self):
 
        print "| test_read_interaction_file | "
 
        self.run.read_transcript_types()
        self.run.read_interaction_file()
 
        #grep O00624 storedInteractions_test.tsv | cut -f1 | cut -f2 -d" " > O00624_11targets
        # grep -f O00624_11targets transcript_biotype.csv 
        # "ENST00000423162","lincRNA","LncRNA"
        # "ENST00000432473","lincRNA","LncRNA"
        # "ENST00000433882","lincRNA","LncRNA"
        # "ENST00000440633","lincRNA","LncRNA"
        # "ENST00000448344","lincRNA","LncRNA"
        # "ENST00000513379","lincRNA","LncRNA"
        # "ENST00000532413","lincRNA","LncRNA"
        # "ENST00000554252","lincRNA","LncRNA"
        # "ENST00000554711","lincRNA","LncRNA"
        # "ENST00000568332","lincRNA","LncRNA"
        # "ENST00000609581","lincRNA","LncRNA"
        
        with open( self.run.outputFile, "r") as inFile:
            header = inFile.readline()
            for line in inFile:
                spl = line.split()
                self.assertTrue( spl[0] == "O00624")
                self.assertTrue( spl[1] == "0")
                self.assertTrue( spl[2] == "11")
                self.assertTrue( spl[4] == "inf")
                break

    # #
    def test_read_sorted_interaction_file(self):
 
        print "| test_sorted_read_interaction_file | "

        # input file needs to be sorted and number of items against MRNAs needs to be the same as against LNCRNAs
        self.interactionsFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/mrna_vs_lncrna/test/proteome_vs_lncrna_mrna_interactions.out_sorted_head"

        self.run = ProteinTargetRatio( self.interactionsFile, self.transcriptTypesFile, self.outputFile)
 
        self.run.read_transcript_types()
 
        self.run.read_sorted_interaction_file()
        
        
        

    def test_fisher_exact_test(self):
        
        print " | test_fisher_exact_test | "
        
        matrix = [ [1,10], [1,10]]
        
        oddsRatio, pvalue = self.run.fisher_exact_test( matrix)

        self.assertTrue( oddsRatio == 1.0)

        matrix = [ [2,10], [1,10]]
        
        oddsRatio, pvalue = self.run.fisher_exact_test( matrix)
        
        self.assertTrue( oddsRatio == 2.0)

        matrix = [ [1,10], [2,10]]
        
        oddsRatio, pvalue = self.run.fisher_exact_test( matrix)
        
        self.assertTrue( oddsRatio == 0.5)
        self.assertTrue( pvalue == 1)
        
        matrix = [ [50,100], [10,100]]
         
        oddsRatio, pvalue = self.run.fisher_exact_test( matrix)
         
        self.assertTrue( pvalue < 0.05)


    def test_t_test(self):
        
        print " | test_t_test | "

        sample1 = [-23,34,0.45,23,-1,3,185,-45.3,-23,34,0.45,23,-1,3,185,-45.3]
        sample2 = [67,45,35,10,56,135,185,45.3,67,45,35,10,56,135,185,45.3]

        # In R
        #             Welch Two Sample t-test
        # 
        # data:  sample1 and sample2
        # t = -2.2807, df = 28.978, p-value = 0.0301

        statistic, pvalue = self.run.t_test(sample1, sample2)
        
        self.assertTrue( statistic - -2.28 < 0.01)
        self.assertTrue( pvalue - 0.0301 < 0.01)


    def test_correct_pvalues(self):

        print " | test_correct_pvalues | "

        dataStore = [["u1"],["u2"],["u3"],["u4"],["u2"],["u3"],["u4"],["u2"],["u3"],["u4"],["u2"],["u3"],["u4"]]
        pvalues = [0.1,0.002,0.001,0.004,0.002,0.001,0.004,0.002,0.01,0.04,0.02,0.1,0.04]

        processedData = self.run.correct_pvalues( len( dataStore), pvalues, dataStore)

        self.assertTrue( processedData[0][-1] == "0")
        self.assertTrue( processedData[1][-1] == "1")
        self.assertTrue( float( processedData[1][-2]) - 7.4e-03 < 0.01)
     
#     # #
#     # Runs after each test
#     def tearDown(self):
#                                 
#         # Wipe output folder
#         cmd = "rm %s/*" % self.outputFolder
#         os.system(cmd)
            
       


