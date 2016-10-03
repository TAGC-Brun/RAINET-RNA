import argparse
import glob
import numpy as np
import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

from fr.tagc.rainet.core.execution.analysis.NetworkScoreAnalysis.NetworkScoreAnalysis import NetworkScoreAnalysis

#===============================================================================
# Started 03-Oct-2016 
# Diogo Ribeiro
DESC_COMMENT = "Combine results from NetworkScoreAnalysis ran for several top numbers."
SCRIPT_NAME = "combine_network_analysis_results.py"
#===============================================================================

# Folder name
FOLDER_NAME = "topPartners"

# Output header
TOP_PARTNERS_HEAD = "TopPartners"
LCN_REAL_HEAD = "LCN_real"
LCN_RANDOM_HEAD = "LCN_random"
SP_REAL_HEAD = "SP_real"
SP_RANDOM_HEAD = "SP_random"

# # column index (0-based) of input files
# TOP_PARTNERS_COL = 0
# LCN_REAL_COL = 1
# LCN_RANDOM_COL = 2
# SP_REAL_COL = 4
# SP_RANDOM_COL = 5


# Read input file and return averages by column
def read_input_file( input_file):
    
    # E.g. 
    #     transcriptID    LCneighbours    LCneighboursRandom      LCneighboursPval        ShortestPath    ShortestPathRandom      ShortestPathPval
    #     ENST00000477643 9       4.39    0.0e+00 3.32    3.81    0.0e+00
    #     ENST00000548215 0       1.98    1.0e+00 4.86    4.27    9.8e-01

    with open( input_file) as inFile:
        
        table = pd.read_table( inFile, header = 0, sep = "\t", skip_blank_lines = True)
        
        lcnRealMean = "%.2f" % np.mean( table["LCneighbours"])
        lcnRandomMean = "%.2f" % np.mean( table["LCneighboursRandom"])
        spRealMean = "%.2f" % np.mean( table["ShortestPath"])
        spRandomMean = "%.2f" % np.mean( table["ShortestPathRandom"])


    return lcnRealMean, lcnRandomMean, spRealMean, spRandomMean
        

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
        parser.add_argument('inFolder', metavar='inFolder', type=str,
                             help='Folder with several topPartners* folders to be analysed.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        filesToProcess = glob.glob( args.inFolder + "/" + FOLDER_NAME + "*/" + NetworkScoreAnalysis.REPORT_METRICS_OUTPUT)
        
        #===============================================================================
        # Write output file
        #===============================================================================
        
        outFile = open("combined_results.tsv", "w")
         
        outFile.write("%s\n" % ( "\t".join( [TOP_PARTNERS_HEAD, LCN_REAL_HEAD, LCN_RANDOM_HEAD, SP_REAL_HEAD, SP_RANDOM_HEAD])) )
        
        resultsDict = {}
        
        for inputFile in filesToProcess:
            #e.g. /TAGC/rainetDatabase/results/networkAnalysis/NetworkScoreAnalysis/100tx_produce_plots/topPartners20/metrics_per_rna.tsv
            topPartnersValue = int( inputFile.split("/")[-2].split(FOLDER_NAME)[1])

            resultsDict[ topPartnersValue] = read_input_file( inputFile) #[ lcnReal, lcnRandom, spReal, spRandom]
            
        for topValue in sorted( resultsDict):
            outFile.write( "%s\t%s\n" % ( topValue, "\t".join( resultsDict[ topValue]) ) )

        outFile.close()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

