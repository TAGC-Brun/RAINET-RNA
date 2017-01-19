
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

from fr.tagc.rainet.core.execution.analysis.EnrichmentAnalysis.LncRNAGroupOddsRatio import LncRNAGroupOddsRatio

# #
# Unittesting the ReadCatrapid script on a specific validated dataset. 
#
class LncRNAGroupOddsRatioUnittest(unittest.TestCase):

    # Constants with default paramters        
        
    # #
    # Runs before each test
    # name of this function needs forcely to be 'setUp'
    def setUp(self):

        # Set the options

        self.annotationFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/lncRNA_groups/transcript_to_gene/lncRNA_group_annotation.txt"
        self.externalFiles = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/test/externalFiles.txt"
        self.backgroundList = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/list_interacting_lncRNAs.txt"
        self.outputFolder = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/analysis/test_output/outfile.tsv"
        self.useGenes = 1
        self.rainetDB = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_backup/RNA/rainet2016-12-28.human_lncRNA_cutoff50/rainet2016-12-28.human_lncRNA_cutoff50.sqlite"        
        
        
        self.run = LncRNAGroupOddsRatio( self.annotationFile, self.externalFiles, self.backgroundList, self.outputFolder, self.useGenes, self.rainetDB)
        

    @unittest.skip("skip cumbersome test")
    def test_read_rainet_db(self):

        print "| test_read_rainet_db | "

        self.run.read_rainet_db()
        
        
        self.assertTrue( len( self.run.rnaCrossReference) == 216133) #number genes = 65988
        

    def test_read_annotation_file(self):
 
        print "| test_read_annotation_file | "
 
        self.run.read_annotation_file()
         
         
        self.assertTrue( len( self.run.groupTranscripts) == 4) 
        self.assertTrue( len( self.run.transcriptAnnotation) ==  12233 )

        #TODO: same test but converting transcript ID to gene ID!

        #diogo@diogo-tower:~/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/Cutoff50/merge/enrichment_lnc2cancer$ cut -f1 /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/lncRNA_groups/transcript_to_gene/lncRNA_group_annotation.txt | sort -u | wc -l
        #12233


    def test_read_annotation_file_transcript_id(self):
 
        print "| test_read_annotation_file_transcript_id | "

        self.run.annotationFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/test/test_annotation_file_transcriptIDs.txt"

        # this should convert transcript ID to gene ID
        self.run.read_annotation_file()
        
        self.assertTrue( len( self.run.groupTranscripts) == 4) 
        self.assertTrue( len( self.run.transcriptAnnotation) ==  3 )

        # same test without conversion
        self.run.useGenes = 0
        self.run.read_annotation_file()

        self.assertTrue( len( self.run.transcriptAnnotation) ==  6 )


    def test_read_external_files(self):
 
        print "| test_read_external_files | "

        self.run.read_external_files()
                
        self.assertTrue( len( self.run.externalLists) == 3) 
        
        # wc -l /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/lncRNA_groups/transcript_to_gene/mukherjee2016_genes.txt
        # removing header and empty line
        self.assertTrue( len( self.run.externalLists["mukherjee2016_genes"]) == 1972) 

        # converted from transcript to gene
        self.assertTrue( len( self.run.externalLists["lnc2cancer_ensemblID"]) == 381) 
        


    def test_read_background_list(self):
 
        print "| test_read_background_list | "

        self.run.read_annotation_file()
        self.run.read_external_files()

        self.run.read_background_list()
        
        # wc -l      12233 list_non_redundant_genes.txt
        self.assertTrue( len( self.run.backgroundTranscripts) == 12233) 

        #  comparelists list_non_redundant_genes.txt mukherjee2016_genes.txt all_combinations 1
        # file1 ids: 12233
        # file2 ids: 1974
        # Intersection: 1129
        self.assertTrue( len( self.run.filteredExternalLists[ "mukherjee2016_genes"]) == 1129) 
        

    def test_calculate_odds_ratio(self):
 
        print "| test_calculate_odds_ratio | "

        self.run.read_annotation_file()
        self.run.read_external_files()
        self.run.read_background_list()

        self.run.calculate_odds_ratio()

        # TODO: test


    def test_fisher_exact_test(self):

        print "| test_fisher_exact_test | "
        
        oddsRatio, pvalue = self.run.fisher_exact_test( [[8, 2], [1, 5] ])

        self.assertTrue( oddsRatio == 20.0)
        self.assertTrue( abs( pvalue - 0.024) <= 0.01)


        
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
            
       


