#=======================================================================
# CONSTANTS FOR RAINET UNITTESTS
#=======================================================================

DB_PATH = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/rainet_testing_DB.sqlite"

#=======================================================================
# InsertionUnittest
#
# Note that values were always obtained externally to the database, 
# e.g. grep of raw input files, looking at the ".dia" Database model etc
# Test input files: /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data
#=======================================================================

TOTAL_NUMBER_TABLES = 31
TOTAL_COLUMNS_IN_RNA_TABLE = 16
TOTAL_COLUMNS_IN_RNAXREF_TABLE = 3
EXAMPLE_MRNA = "ENST00000379749"
EXAMPLE_MRNA_TABLE_COLUMNS = {
                    "transcriptID": "ENST00000379749",
                    "geneID": "ENSG00000133110",
                    "peptideID": "ENSP00000369073",
                    "transcriptBiotype": "protein_coding",
                    "transcriptLength": 3210,
                    "transcriptSource": "havana",
                    "transcriptStatus": "KNOWN",
                    "transcriptTsl": "tsl5 (assigned to previous version 7)",
                    "transcriptGencodeBasic": 1,
                    "transcriptStart": 37562585,
                    "transcriptEnd": 37598761,
                    "transcriptStrand": -1,
                    "chromosomeName": "13",
                    "percentageGCContent": 32.94,
                    "description" : "periostin, osteoblast specific factor [Source:HGNC Symbol;Acc:HGNC:16953]",
                    "type": "MRNA"
                    }
EXAMPLE_MRNA_XREFS = {
                     "hgnc_id|HGNC:16953",
                     "hpa|HPA012306",
                     "ucsc|uc058wmq.1",
                     "embl|AL138679",
                     "embl|AL646087",
                     "uniprot_sptrembl|B1ALD9",
                     "refseq_mrna_predicted|XM_005266231"
                    }
EXAMPLE_MRNA_PEPTIDE = ["ENST00000379749","B1ALD9"]
EXAMPLE_PRI_RNA = "ENST00000411230"
EXAMPLE_PRI_RNA_INTERACTIONS = {
                                "ENSP00000410728" : -7.61,
                                "ENSP00000450122" : 4.02,
                                "ENSP00000431009" : 8.68
                                }
EXAMPLE_PRI_RNA_PROTEIN_UNIPROT = {"ENSP00000410728" : "Q5RI18"} 
NUMBER_TISSUES = 49

#=======================================================================

#=======================================================================
# AnalysisStrategyUnittest
#
#=======================================================================

