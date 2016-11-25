import sys
import os
import argparse
import glob

# import numpy as np
import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil

from constants import *

#===============================================================================
# Started 24-Nov-2016 
# Diogo Ribeiro
DESC_COMMENT = "Wrapper script for ENCODE eCLIP dataset. To filter, merge replicates, map against gene models, produce interactions file."
SCRIPT_NAME = "process_all_eclip_files.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read metadata file, have group replicates, match identifiers
# 2) Read bed files, apply filtering
# 3) Merge replicate bed files
# 4) Map peaks to transcripts
# 5) Produce interaction file
#===============================================================================

#===============================================================================
# Processing notes:
# 1)
# 2)
#===============================================================================


FILTER_ECLIP_FILES_SCRIPT_PATH = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/knownScaffoldExamples/eCLIP/filter_eclip_files.py"


# #
# Read, filter and write filtered bed file
def read_bed_folder( bedFolder, minPval, minFC):

    filesToProcess = glob.glob( bedFolder + "/*.bed" )

    print "Processing %s files.." % len( filesToProcess)
    
    
    for fi in filesToProcess:
        cmd = "workon rainet; python %s %s --minPval %s --minFC %s" % ( FILTER_ECLIP_FILES_SCRIPT_PATH, fi, minPval, minFC)

        returnCode = SubprocessUtil.run_command( cmd)

        #TODO: here
        print fi


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
        parser.add_argument('bedFolder', metavar='bedFolder', type=str,
                             help='Input folder with all bed files to be processed.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
        parser.add_argument('metadataFile', metavar='metadataFile', type=str,
                             help='Metadata file provided when downloding the bed files.')
        parser.add_argument('bedModels', metavar='bedModels', type=str,
                             help='Bed file with gene models for identifying overlapping transcripts.')
        parser.add_argument('--minPval', metavar='minPval', type=float, default = MIN_PVAL_DEFAULT,
                             help='Peaks below given value will be excluded. Note that provided pvalue is positive log10, the higher the value, the more significant it is. (Default = -1, i.e. "OFF").')
        parser.add_argument('--minFC', metavar='minPval', type=float, default = MIN_FC_DEFAULT,
                             help='Peaks below given value will be excluded. Note that provided fold-change enrichment is positive log2, the higher the value, the higher is the fold change. (Default = -1, i.e. "OFF").')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        if not os.path.exists( args.outputFolder):
            os.mkdir( args.outputFolder)

        Timer.get_instance().step( "Read bed folder..")            
        read_bed_folder( args.bedFolder, args.minPval, args.minFC)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

