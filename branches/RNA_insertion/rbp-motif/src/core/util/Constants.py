
import os



PATH_HOME = os.path.expanduser('~')
PATH_LOG = os.path.expanduser('~') + '/workspace/tagc-rainet/rbp-motif/output/log_file/RBP_motif.log'
PATH_ANCHOR = os.environ.get('ANCHOR_PATH')
PATH_IUPRED = os.environ.get('IUPred_PATH')

#==============================================================================
# Constant file name for InfoDataset
#==============================================================================

FILE_COMMON = 'Common_items.txt'
FILE_DIFF = 'Diff_items.txt'

#===============================================================================
# Constants for Logging 
#===============================================================================
MODE_DEBUG = "debug"
MODE_INFO = "info"
MODE_WARNING = "warning"
MODE_ERROR = "error"
MODE_CRITICAL = "critical"
MODE_FATAL = "fatal"
VERBOSITY_LEVELS = [ MODE_DEBUG, MODE_INFO, MODE_WARNING, MODE_ERROR, MODE_CRITICAL, MODE_FATAL]

#LOG_MODE_DEBUG = ['d', 'debug']
# LOG_MODE_INFO = ['i', 'info']
# LOG_MODE_WARNING = ['w', 'warning']
# LOG_MODE_ERROR = ['e', 'error']
# LOG_MODE_CRITICAL = ['c', 'critical']


LOG_APPEND = 'a'
LOG_NO_APPEND = 'w'


#===============================================================================
# Extensions files tag
#===============================================================================
EXTENSION_ODS_TAG = 'ODS'
EXTENSION_XLS_TAG = 'XLS'
EXTENSION_XLSX_TAG = 'XLSX'
EXTENSION_TXT_TAG = 'TXT'
EXTENSION_CSV_TAG = 'CSV'
EXTENSION_LIST = [EXTENSION_ODS_TAG, EXTENSION_XLS_TAG,
                  EXTENSION_XLSX_TAG, EXTENSION_TXT_TAG, EXTENSION_CSV_TAG]

#===============================================================================
# Constants on data insertion status 
#===============================================================================
STATUS_OK = "OK"
STATUS_WARNING = "WARNING"
STATUS_RAINET_ERROR = "RAINET ERROR"
STATUS_ERROR = "ERROR" 

#===============================================================================
# Constants for connection to Ensembl
#===============================================================================


SERVER_RELEASE_75 = "http://feb2014.archive.ensembl.org"
END_URL_ENSEMBL = ";output=fasta;param=peptide;_format=Text"
SECTION_ENSEMBL = "/Homo_sapiens/Export/Output/Gene?"
TYPE_QUERY_ENSEMBL = {'ENSG':'g=', 'ENST':'t='}


#===============================================================================
# Constants for connection to Uniprot
#===============================================================================


SERVER_UNIPROT = 'http://www.uniprot.org/uniprot/'
SEQ = '.fasta'
INFO = '.txt'
END_URL_UNIPROT = [INFO, SEQ]


# ==============================================================================
# Web connection constant
# ==============================================================================

WEB_CONNECTION_CONSTANT = 5



# ==============================================================================
# DisoRDPbind constant
# ==============================================================================

BINDING_PARTNER = {1:'RNA-binding residues', 2:'DNA-binding residues', 3:'protein-binding residues'}


# ===============================================================================
# ClassificationNewDataset
# ===============================================================================

DOMAIN_CLASSICAL = 'Classical'
DOMAIN_NONCLASSICAL = 'Non-classical'
DOMAIN_UNKNOWN = 'Unclassified'
DOMAIN_NONE = 'No-domains'
DOMAIN_OTHER = 'Other-Domain'
DOMAIN_COMMA = ','



#=================================================================================
# OverlapMotifsRegion
#================================================================================

ELM_FILE = 'FMotifs_'
EXTENSION_TXT_FILE = '.txt'
IUPRED_PRED_FILE = 'IUPred_'
ANCHOR_PRED_FILE = "Pred_"
THRESHOLD_LOW = 0.4
THRESHOLD_HIGH = 0.5
OVERLAP_YES = 'Y'
OVERLAP_NO = 'N'
KEY_IUPRED = "iupred table"
KEY_ANCHOR = "anchor table"
KEY_DISORDP = "DisoRDPbind table"

#=================================================================================
# GlobalOVerlapMotifsRegion
#=================================================================================

IUPRED_TOOL = 'IUPRED'
ANCHOR_TOOL = 'ANCHOR'
DISORDP_TOOL = 'DISORDP'
DICT_FILE_OUTPUT  = {'IUPRED': ('SLimIUPred_Overlap_','SLimDomainIUPred_Overlap_'), 'ANCHOR': ('SLimANCHOR_Overlap_','SLimDomainANCHOR_Overlap_'),
                     'DISORDP':('SLimDisoRDP_Overlap_', 'SLimDomainDisoRDP_Overlap')}

MOTIFS_THRESHOLD_1 = 0.4
MOTIFS_THRESHOLD_2 = 0.5

SLIM_NAME = 'SLiM NAME'
SLIM_REGION_COUNTER = 'SLiM Region Number'
START_SLIM_REGION = 'Start'
END_SLIM_REGION = 'End'
TOOLNAME = 'Prediction Tool'
TOOLNAME_REGION_COUNTER = 'Tool Region Number'
START_TOOL_REGION = 'Start'
END_TOOL_REGION = 'End'
OVERLAP_OUTCOME = 'Overlap Outcome'
OVERLAP_LENGTH = 'Overlap Length'
 

