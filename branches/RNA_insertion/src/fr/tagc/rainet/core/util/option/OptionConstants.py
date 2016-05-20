
import os

from fr.tagc.rainet.core.util import Constants
from fr.tagc.rainet.core.data import DataConstants

#===============================================================================
# The list of available strageties
#===============================================================================
STRATEGIES_LIST = [ "Insertion", "InteractiveQuery", "Analysis", "DatabaseCheck"]

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
OPTION_OUTPUT_FOLDER = "Output folder"

#===============================================================================
# Constants for analysis strategy
#===============================================================================

DEFAULT_BIOTYPE = "RNA"
DEFAULT_INTERACTION_SCORE = "OFF"
DEFAULT_GENCODE = 0
DEFAULT_EXPRESSION_VALUE_CUTOFF = 0
DEFAULT_OUTPUT_FOLDER = os.getcwd()

RNA_BIOTYPES = DataConstants.RNA_MRNA_BIOTYPE[:] + DataConstants.RNA_LNCRNA_BIOTYPE[:]
DEFAULT_RNA_BIOTYPES = RNA_BIOTYPES[:]

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
                    [ "-e", "--expressionValueCutoff", "store", "float", OPTION_EXPRESSION_VALUE_CUTOFF, DEFAULT_EXPRESSION_VALUE_CUTOFF, "Protein-RNA interactions where one of its components has expression below the given value will be excluded. Default: OFF"]
                ],
                "DatabaseCheck" : [
                    [ "-d", "--databasePath", "store", "string", OPTION_DB_NAME, None, "The path to the SQL database to use/create."],
                    [ "-v", "--verbose", "store", "string", OPTION_VERBOSITY, Constants.MODE_INFO, "The level of verbosity. Must be one of : " + str( Constants.VERBOSITY_LEVELS)]
                ]
               }
