import sys
import os
import argparse

# import numpy as np
# import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

# from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
# from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
# from fr.tagc.rainet.core.util.data.DataManager import DataManager

# from fr.tagc.rainet.core.data.Protein import Protein


#===============================================================================
# Started 25-June-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to ..."
SCRIPT_NAME = "name_of_script.py"
#===============================================================================

#===============================================================================
# General plan:
# 1)
# 2)
#===============================================================================

#===============================================================================
# Processing notes:
# 1)
# 2)
#===============================================================================

# #
# Write line about what the function does
def function_one( arg1):

#     SQLManager.get_instance().set_DBpath( DBPATH)
#     sql_session = SQLManager.get_instance().get_session()

    pass        


if __name__ == "__main__":

    try:
    
        # Start chrono
        Timer.get_instance().start_chrono()
        print "STARTING " + SCRIPT_NAME
        
        #===============================================================================
        # Get input arguments
        #===============================================================================
        parser = argparse.ArgumentParser(description= DESC_COMMENT) 
    
        # positional args
        parser.add_argument('sysArgFile', metavar='sysArgFile', type=str,
                             help='Description.')
        parser.add_argument('--sysArgOptional', metavar='sysArgOptional', type=int, default = 1,
                             help='Description (Default = 1).')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        # run function to..
        Timer.get_instance().step( "Read bogus file..")            
        function_one( args.enrichmentPerRNAFile)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

