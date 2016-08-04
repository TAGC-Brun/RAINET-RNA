
import os

from fr.tagc.rainet.core.util import Constants
from fr.tagc.rainet.core.data import DataConstants

#===============================================================================
# The list of available strageties
#===============================================================================
STRATEGIES_LIST = [ "Insertion", "InteractiveQuery", "Analysis", "DatabaseCheck", "EnrichmentAnalysis"]

#===============================================================================
# The tags used to define the options
#===============================================================================
SHORT_TAG = "short"
LONG_TAG = "long"
ACTION_TAG = "action"
TYPE_TAG = "type"
DEST_TAG = "dest"
DEFAULT_TAG = "default"
HELP_TAG = "help"

#===============================================================================
# The name of the various options
#===============================================================================
# Common options
OPTION_DB_NAME = "Database path"
OPTION_SPECIES = "Species"
OPTION_VERBOSITY = "Verbosity"
OPTION_OUTPUT_FOLDER = "Output folder"
# Insertion strategy options
OPTION_INSERTION_PROPERTIES_PATH = "insertionPropertiesPath"
OPTION_INSERTION_FORCE_OVERRIDE = "insertionForceOverride"
# InteractiveQuery Strategy options
OPTION_QUERY_FILE = "query_file"
# Analysis Strategy options
OPTION_MINIMUM_INTERACTION_SCORE = "Minimum interaction score"
OPTION_RNA_BIOTYPES = "RNA biotypes filtering"
OPTION_GENCODE = "Gencode basic filtering"
OPTION_EXPRESSION_VALUE_CUTOFF = "Minimum expression value"
OPTION_EXPRESSION_TISSUE_CUTOFF = "Minimum number tissues co-present"
# Enrichment analysis options
OPTION_ANNOTATION_TABLE = "Protein annotation"
OPTION_MINIMUM_PROTEIN_ANNOTATION = "Minimum proteins in annotation"
OPTION_MINIMUM_PROTEIN_INTERACTION = "Minimum proteins in annotation with interaction"
OPTION_NUMBER_RANDOMIZATIONS = "Number randomizations"
OPTION_LOW_MEMORY = "Low memory"

#===============================================================================
# Constants for default values
#===============================================================================

# Analysis strategy
DEFAULT_BIOTYPE = "RNA"
DEFAULT_INTERACTION_SCORE = "OFF"
DEFAULT_GENCODE = 0
DEFAULT_EXPRESSION_VALUE_CUTOFF = "OFF"
DEFAULT_EXPRESSION_TISSUE_CUTOFF = 1.0
DEFAULT_OUTPUT_FOLDER = os.getcwd()
DEFAULT_LOW_MEMORY = 0

RNA_BIOTYPES = DataConstants.RNA_ALL_BIOTYPE # DataConstants.RNA_MRNA_BIOTYPE[:] + DataConstants.RNA_LNCRNA_BIOTYPE[:]
DEFAULT_RNA_BIOTYPES = RNA_BIOTYPES[:]

# Enrichment strategy
DEFAULT_MINIMUM_PROTEIN_ANNOTATION = 5
DEFAULT_MINIMUM_PROTEIN_INTERACTION = 2
DEFAULT_NUMBER_RANDOMIZATIONS = 100

#===============================================================================
# The definition of the options
#===============================================================================
OPTION_LIST = {  "Insertion": [
                    [ "-d", "--databasePath", "store", "string", OPTION_DB_NAME, None, "The path to the SQL database to use/create."],
                    [ "-s", "--species", "store", "string", OPTION_SPECIES, None, "The species used in the database."],
                    [ "-v", "--verbose", "store", "string", OPTION_VERBOSITY, Constants.MODE_INFO, "The level of verbosity. Must be one of : " + str( Constants.VERBOSITY_LEVELS)],
                    [ "-i", "--insertionPropertiesPath", "store", "string", OPTION_INSERTION_PROPERTIES_PATH, None, "The path to the properties file containing the list of files to use for DB insertion."],
                    [ "-f", "--forceOverride", "store_true", None, OPTION_INSERTION_FORCE_OVERRIDE, None, "Indicates if the whole database must be dropped and re-inserted or not."],
                ], 
                 "InteractiveQuery":[
                    [ "-d", "--databasePath", "store", "string", OPTION_DB_NAME, None, "The path to the SQL database to use/create."],
                    [ "-s", "--species", "store", "string", OPTION_SPECIES, None, "The species used in the database."],
                    [ "-q", "--queryFilePath", "store", "string", OPTION_QUERY_FILE, None, "The path to the file containing the query definition."],                 
                ],               
                "Analysis" : [
                    [ "-d", "--databasePath", "store", "string", OPTION_DB_NAME, None, "The path to the SQL database to use/create."],
                    [ "-s", "--species", "store", "string", OPTION_SPECIES, None, "The species used in the database."],
                    [ "-v", "--verbose", "store", "string", OPTION_VERBOSITY, Constants.MODE_INFO, "The level of verbosity. Must be one of : " + str( Constants.VERBOSITY_LEVELS)],
                    [ "-o", "--outputFolder", "store", "string", OPTION_OUTPUT_FOLDER, DEFAULT_OUTPUT_FOLDER, "Folder where output files will be written."],
                    [ "-m", "--minimumInteractionScore", "store", "string", OPTION_MINIMUM_INTERACTION_SCORE, DEFAULT_INTERACTION_SCORE, "Protein-RNA interactions with interaction score below given value will be excluded from analysis. Default: OFF"],
                    [ "-b", "--RNABiotype", "store", "string", OPTION_RNA_BIOTYPES, DEFAULT_RNA_BIOTYPES, "Comma-separated list of RNA biotypes to be INCLUDED in analysis. Default: all biotypes are considered. Can one or several of :"+ str(RNA_BIOTYPES)],
                    [ "-g", "--gencodeBasicOnly", "store", "int", OPTION_GENCODE, DEFAULT_GENCODE, "If 1, include in analysis ONLY transcripts tagged as present in 'GENCODE basic'. Default: 0 (i.e. all RNAs are considered)."],
                    [ "-e", "--expressionValueCutoff", "store", "string", OPTION_EXPRESSION_VALUE_CUTOFF, DEFAULT_EXPRESSION_VALUE_CUTOFF, "Protein-RNA interactions where one of its components has expression below the given value will be excluded. Default: 0"],
                    [ "-t", "--expressionTissueCutoff", "store", "float", OPTION_EXPRESSION_TISSUE_CUTOFF, DEFAULT_EXPRESSION_TISSUE_CUTOFF, "Protein-RNA interactions between pairs co-present in less that this numbers of tissues will be excluded. Only active if expressionValueCutoff is also active. Default: 1"],
                    [ "-l", "--lowMemory", "store", "int", OPTION_LOW_MEMORY, DEFAULT_LOW_MEMORY, "If 1, use less memory but do not produce a full report. Used for producing file with expression filtering. Default: 0"]
                ],
                "EnrichmentAnalysis" : [
                    [ "-d", "--databasePath", "store", "string", OPTION_DB_NAME, None, "The path to the SQL database to use/create."],
                    [ "-s", "--species", "store", "string", OPTION_SPECIES, None, "The species used in the database."],
                    [ "-v", "--verbose", "store", "string", OPTION_VERBOSITY, Constants.MODE_INFO, "The level of verbosity. Must be one of : " + str( Constants.VERBOSITY_LEVELS)],
                    [ "-o", "--outputFolder", "store", "string", OPTION_OUTPUT_FOLDER, DEFAULT_OUTPUT_FOLDER, "Folder where output files will be written."],
                    [ "-a", "--annotationTable", "store", "string", OPTION_ANNOTATION_TABLE, None, "Protein annotation table to be used for analysis. Must be one of : " + str( Constants.ANNOTATION_TABLES) ],
                    [ "-m", "--minimumProteinAnnotation", "store", "int", OPTION_MINIMUM_PROTEIN_ANNOTATION, DEFAULT_MINIMUM_PROTEIN_ANNOTATION, "Minimum number of proteins with annotation for enrichment test to be performed." ],
                    [ "-i", "--minimumProteinInteraction", "store", "int", OPTION_MINIMUM_PROTEIN_INTERACTION, DEFAULT_MINIMUM_PROTEIN_INTERACTION, "Minimum number of proteins in a given annotation with positive interactions for enrichment test to be performed." ],
                    [ "-r", "--numberRandomizations", "store", "int", OPTION_NUMBER_RANDOMIZATIONS, DEFAULT_NUMBER_RANDOMIZATIONS, "Number of randomizations to be performed for the control experiment." ]
                ],
                "DatabaseCheck" : [
                    [ "-d", "--databasePath", "store", "string", OPTION_DB_NAME, None, "The path to the SQL database to use/create."],
                    [ "-v", "--verbose", "store", "string", OPTION_VERBOSITY, Constants.MODE_INFO, "The level of verbosity. Must be one of : " + str( Constants.VERBOSITY_LEVELS)]
                ]
               }
