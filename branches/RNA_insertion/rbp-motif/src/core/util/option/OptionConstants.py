#from core.util import Constants


from core.util import Constants
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
# The name of the various options DA MODIFICIARE
#===============================================================================
OPTION_INPUT_PATH = "--inputPath"
OPTION_OUTPUT_PATH = "--outputPath"
OPTION_FILENAME = "--filename"
OPTION_MORE_FILENAME = "--MoreFilename"
OPTION_NUMBER_INDEX = "--index"
OPTION_NUMBER_LENGTH = "--length"
OPTION_MAKEDATASETRBP_PROPERTIES_PATH = "--MakeDatasetRbpPath"
OPTION_DOWNLOADENSEMBLSEQ_PROPERTY_PATH = "--DownloadEnsemblSeq"
OPTION_DISORDER_ANALYSIS_PATH = "--DisorderAnalysis"
OPTION_MOTIFS_ANALYSIS_PROPERTY_PATH = "--MotifsAnalysis"
OPTION_VERBOSITY = "--verbose"



#===============================================================================
# The definition of the options
#===============================================================================
OPTION_LIST = [ [ "-i", "--inputPath", "store", "string", OPTION_INPUT_PATH, None, " Directory path that will be contain the input file (without the HOME path)", 1],
                [ "-o", "--outputPath", "store", "string", OPTION_OUTPUT_PATH, None, " Directory path of output file (without the HOME path)", 1],
                [ "-f", "--filename", "store", "string", OPTION_FILENAME, None, " File name (without the whole path)", 1],
                ["-v", "--verbose", "store", "string", OPTION_VERBOSITY, Constants.MODE_INFO, "The level of verbosity. Must be one of : " + str( Constants.VERBOSITY_LEVELS), 1],
                [ "-k", "--MoreFilename", "store", "string", OPTION_MORE_FILENAME, None, " You must provide two file name (without the whole path)", 2],
                [ "-n", "--index", "store", "int", OPTION_NUMBER_INDEX, None, " Integer Numbers representing a index of table", 2],
                [ "-l", "--length", "store", "int", OPTION_NUMBER_LENGTH, None, " Integer Numbers representing a length of list", 2],
                [ "-p", "--MakeDatasetRbpPath", "store", "string", OPTION_MAKEDATASETRBP_PROPERTIES_PATH, None, " Property File Path", 1],
                [ "-w", "--DownloadEnsemblSeq", "store", "string" , OPTION_DOWNLOADENSEMBLSEQ_PROPERTY_PATH, None , " Property File Path", 1],
                [ "-d", "--DisorderAnalysis", "store", "string" , OPTION_DISORDER_ANALYSIS_PATH, None , " Property File Path", 1],
                [ "-m", "--MotifsAnalysis", "store", "string", OPTION_MOTIFS_ANALYSIS_PROPERTY_PATH,None, "Property file path", 1]
                ]
