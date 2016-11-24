import sys
import os
import argparse

# import numpy as np
import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

# from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil

#===============================================================================
# Started 24-Nov-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to filter ENCODE eCLIP bed files based on pvalue or fold-enrichment."
SCRIPT_NAME = "filter_eclip_files.py"
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

# 0-based bed columns
PVALUE_COLUMN = 7
FOLD_CHANGE_COLUMN = 6

# previous file format (E.g. 293T cell line)
# PVALUE_COLUMN = 3
# FOLD_CHANGE_COLUMN = 4

# #
# Read, filter and write filtered bed file
def read_bed_file( bedFile, minPval, minFC):

    # Example format
    # chr7    128862938       128863050       FTO_K562_rep01  1000    +       3.0648192732971 26.9230684254443        -1      -1
    # chr7    140696769       140696843       FTO_K562_rep01  1000    +       4.01319974942554        25.348382825905 -1      -1
    # whereas the (0-based) column 6 is the pvalue and column 7 the fold enrichment.

    # pvalue is: log10(p-value eCLIP vs SMInput). Calculated by: my $log10pval = $chipval > 0 ? -1 * log($chipval)/log(10) : 400 ;
    # therefore, the higher, the more significant

    # fold-change is: log2(fold-enrichment in eCLIP vs SMInput). Calculated by: my $l2fc = log(($peak_read_counts{$peak}{"expt"}/$mapped_read_count{"expt"}) / ($peak_read_counts{$peak}{"input"}/$mapped_read_count{"input"})) / log(2);
    # therefore, the higher the value, the higher is the fold change

    table = pd.read_table( bedFile, header = None, sep = "\t", skip_blank_lines = True, skiprows = 0)

    print "Initial number of lines in file:", len( table)

    filteredTable = table.copy()

    ### Pvalue filter
    # if different from defualt value
    if minPval != -1:
        filteredTable = filteredTable.loc[table.loc[:, PVALUE_COLUMN] > minPval,:]

    print "After pvalue filter:", len( filteredTable)

    ### Fold-change enrichment filter
    # if different from defualt value
    if minFC != -1:
        filteredTable = filteredTable.loc[table.loc[:, FOLD_CHANGE_COLUMN] > minFC,:]

    print "After pvalue and Fold-change enrichment filter:", len( filteredTable)

    ### write filtered table to output file
    pd.DataFrame.to_csv( filteredTable, bedFile + "_filtered", sep = "\t", header = False, index = False)
    #Note: float precision is changed when using pandas

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
        parser.add_argument('bedFile', metavar='bedFile', type=str,
                             help='Input bed file to be processed.')
        parser.add_argument('--minPval', metavar='minPval', type=float, default = -1,
                             help='Peaks below given value will be excluded. Note that provided pvalue is positive log10, the higher the value, the more significant it is. (Default = -1, i.e. "OFF").')
        parser.add_argument('--minFC', metavar='minPval', type=float, default = -1,
                             help='Peaks below given value will be excluded. Note that provided fold-change enrichment is positive log2, the higher the value, the higher is the fold change. (Default = -1, i.e. "OFF").')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        Timer.get_instance().step( "Read bed file..")            
        read_bed_file( args.bedFile, args.minPval, args.minFC)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

