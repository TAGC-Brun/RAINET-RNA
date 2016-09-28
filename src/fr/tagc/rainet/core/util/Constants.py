
import os

PATH_LOG = os.path.expanduser('~') + '/.rainet/Rainet.log'


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
# 
# LOG_MODE_DEBUG = ['d', 'debug']
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
# Constants on folder creation
#===============================================================================
REPORT_FOLDER = "Report"

#===============================================================================
# Constants on enrichment analysis
#===============================================================================
ANNOTATION_TABLES = ["NetworkModule", "ReactomePathway", "KEGGPathway", "BioplexCluster", "CorumCluster", "WanCluster", "CustomCluster"]


