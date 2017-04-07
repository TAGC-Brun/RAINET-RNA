
import unittest
import os
import pandas as pd
import glob

from fr.tagc.rainet.core.execution.analysis.EnrichmentAnalysis.CommonLncRNAProteinDisease import CommonLncRNAProteinDisease

# #
# Unittesting 
#
class CommonLncRNAProteinDiseaseUnittest(unittest.TestCase):

    # Constants with default paramters        
        
    # #
    # Runs before each test
    # name of this function needs forcely to be 'setUp'
    def setUp(self):

        print os.getcwd()

        # Set the options
        self.rainetDB = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_backup/RNA/rainet2016-12-28.human_lncRNA_cutoff50/rainet2016-12-28.human_lncRNA_cutoff50.sqlite"
        self.lncRNADiseaseFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_input/EnrichmentAnalysis/CommonLncRNAProteinDisease/enriched_lncRNA_disease_data_lncrnadisease_lnc2cancer.txt"
        self.proteinDiseaseFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_input/EnrichmentAnalysis/CommonLncRNAProteinDisease/protein_disease_description.txt"
        self.enrichmentData = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_input/EnrichmentAnalysis/CommonLncRNAProteinDisease/list_filtered_enrichments_rnas_havugimana2012_bioplex_wan2015_networkModules.txt"
        self.outputFolder = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/enrichmentAnalysis/CommonLncRNAProteinDisease"
        self.minWordSize = 2
        self.blackListedWords = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_input/EnrichmentAnalysis/CommonLncRNAProteinDisease/black_listed_words.txt"
        self.complexDatasets = "WanCluster,BioplexCluster,CustomCluster,NetworkModule"

        # folder containing expected output files
        self.expectedFolder = "test_expected/"
        
        self.run = CommonLncRNAProteinDisease(self.rainetDB, self.lncRNADiseaseFile, self.proteinDiseaseFile, self.enrichmentData, 
                                self.outputFolder, self.minWordSize, self.blackListedWords, self.complexDatasets )
            

    # #
    def test_enrichment_to_protein(self):

        print "| test_enrichment_to_protein | "

        pairsToWrite, enrichmentsDict = self.run.enrichment_to_protein()
        
        # ENST00000320202 vs Bioplex cluster 185 (O95721, Q96NA8, Q9NXR1)
    
        self.assertTrue( len(pairsToWrite) == 267494)
        
        self.assertTrue( "ENST00000320202\tO95721\n" in pairsToWrite)
        
        self.assertTrue( len( enrichmentsDict) == 27090)

        # ENST00000564147    250    27    37    40.2%    0    1.2e-05    1.1e-03    1    NetworkModule

        res = enrichmentsDict["ENST00000564147|250|NetworkModule"].split("\t")
        
        self.assertTrue( len( res[-1].split(",")) == 52, "assert that number of proteins in complex is correct")
        self.assertTrue( int( res[2]) == 27)


    def test_process_lncrna_protein_map(self):
        
        print "| test_process_lncrna_protein_map | "
        
        self.run.read_black_listed_words_file()           

        self.run.read_lncrna_disease_file( )

        self.run.read_rainet_db( )
        
        self.run.read_protein_disease_file( )

        self.run.enrichment_to_protein()

        self.run.process_lncrna_protein_map()

        # before I found 127 word match lines, now I find 165. 
        # Probably because now I also have complex level, even if tx-protein pair is the same, but from different complex, they will now be in a different line

        ################### TODO: test this ######################


    def test_read_rainet_db(self):
        
        print "| test_read_rainet_db | "

        self.run.read_lncrna_disease_file( )
        
        self.run.read_rainet_db( )

        self.assertTrue( len( self.run.transcriptInteractionDict) == 35)
        
        self.assertTrue( len( self.run.transcriptInteractionDict["ENST00000432701"]) == 3214)
        
        self.assertTrue( len( self.run.geneTranscriptDict) == 35)

        self.assertTrue( self.run.geneTranscriptDict["ENST00000582727"] == "ENSG00000232956")

     
#     # #
#     # Runs after each test
#     def tearDown(self):
#                                  
#         # Wipe output folder
#         cmd = "rm %s/*" % self.outputFolder
#         os.system(cmd)
            
       


