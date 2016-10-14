import argparse
import glob
import numpy as np
import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

from statsmodels.stats.multitest import multipletests

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

PVAL_THRESHOLD = 0.05

NUMBER_EXAMPLE_TRANSCRIPTS = 20

# Read input file and return averages by column
def read_input_file( input_file, only_significant):
    
    # E.g. 
    #     transcriptID    LCneighbours    LCneighboursRandom      LCneighboursPval        ShortestPath    ShortestPathRandom      ShortestPathPval
    #     ENST00000477643 9       4.39    0.0e+00 3.32    3.81    0.0e+00
    #     ENST00000548215 0       1.98    1.0e+00 4.86    4.27    9.8e-01

    exampleTranscripts = set()
    boo = 1

    with open( input_file) as inFile:
        
        table = pd.read_table( inFile, header = 0, sep = "\t", skip_blank_lines = True)

        # pick example transcript for the first file read
        if boo:
            exampleTranscripts = set( table["transcriptID"][0:NUMBER_EXAMPLE_TRANSCRIPTS])
            boo = 0

        # Get LCN just for the example transcripts
        # order of transcripts is always kept in input files and thus also in the table
        exampleLCNs = list( table.loc[table["transcriptID"].isin(exampleTranscripts)]["LCneighbours"] )
        exampleLCNs = [ "%.2f" % val for val in exampleLCNs]
       
        # whether to keep only transcripts with significant results compared to random
        if only_significant:
            # filter by LCneighbours p value
            if only_significant == 1:        
                cols = ["LCneighboursPval"]
                table[cols] = table[table[cols] < PVAL_THRESHOLD][cols]
                table = table[table.LCneighboursPval.notnull()]
            # filter by Shortest path p value
            if only_significant == 2:
                cols = ["ShortestPathPval"]
                table[cols] = table[table[cols] < PVAL_THRESHOLD][cols]
                table = table[table.ShortestPathPval.notnull()]

        # count number significant entries for LCN
        newTable = table.copy()
        cols = ["LCneighboursPval"]
        # apply multiple test correction
        newTable[cols] = multiple_test_correction( table.LCneighboursPval)
        newTable[cols] = newTable[newTable[cols] < PVAL_THRESHOLD][cols]
        newTable = newTable[newTable.LCneighboursPval.notnull()]        
        signLCN = str( len(newTable))


        # count number significant entries for SP
        newTable = table.copy()
        cols = ["ShortestPathPval"]
        # apply multiple test correction
        newTable[cols] = multiple_test_correction( table.ShortestPathPval)
        newTable[cols] = newTable[newTable[cols] < PVAL_THRESHOLD][cols]
        newTable = newTable[newTable.ShortestPathPval.notnull()]
        signSP = str( len(newTable))

        lcnRealMean = "%.2f" % np.mean( table["LCneighbours"])
        lcnRandomMean = "%.2f" % np.mean( table["LCneighboursRandom"])
        spRealMean = "%.2f" % np.mean( table["ShortestPath"])
        spRandomMean = "%.2f" % np.mean( table["ShortestPathRandom"])
        nTranscripts = str( len( table) ) #n transcripts used for means etc


        ## write filtered table
        table.to_csv(input_file + "_processed.txt", sep = "\t")
        
    return lcnRealMean, lcnRandomMean, spRealMean, spRandomMean, nTranscripts, signLCN, signSP, exampleTranscripts, exampleLCNs
        

# #
# Run multipletest correction and return pvalues
def multiple_test_correction(pvalues, meth = "fdr_bh"):
    return multipletests(pvalues, method = meth)[1]
        

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
        parser.add_argument('--onlySignificant', metavar='onlySignificant', type=int, default = 0,
                             help='If 1, calculate means only when differences are significant. If 2, use shortestPath instead of LCneighbours for the filter.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        filesToProcess = glob.glob( args.inFolder + "/" + FOLDER_NAME + "*/" + NetworkScoreAnalysis.REPORT_METRICS_OUTPUT)
        
        if len( filesToProcess) == 0:
            raise RainetException( "combine_network_analysis_results : No files to process in %s" % ( args.inFolder ) )
        
        #===============================================================================
        # Write output file
        #===============================================================================
        
        outFile = open("combined_results.tsv", "w")
        outFile2 = open("combined_examples.tsv", "w")
         
        outFile.write("%s\n" % ( "\t".join( [TOP_PARTNERS_HEAD, LCN_REAL_HEAD, LCN_RANDOM_HEAD, SP_REAL_HEAD, SP_RANDOM_HEAD, "n_transcripts", "signLCN", "signSP"])) )
        
        resultsDict = {}
        
        exampleDict = {}
        
        for inputFile in filesToProcess:
            #e.g. /TAGC/rainetDatabase/results/networkAnalysis/NetworkScoreAnalysis/100tx_produce_plots/topPartners20/metrics_per_rna.tsv
            topPartnersValue = int( inputFile.split("/")[-2].split(FOLDER_NAME)[1])

            data = read_input_file( inputFile, args.onlySignificant)
            resultsDict[ topPartnersValue] = data[0:-2] #[ lcnReal, lcnRandom, spReal, spRandom]

            # get the LCN values for example transcripts          
            exampleDict[ topPartnersValue] = data[ -1]
            # get list of example transcripts IDs
            exampleHeader = data[ -2]


        outFile2.write("topPartners\t%s\n" % ( "\t".join( exampleHeader)) )
            
        for topValue in sorted( resultsDict):
            outFile.write( "%s\t%s\n" % ( topValue, "\t".join( resultsDict[ topValue]) ) )
            outFile2.write( "%s\t%s\n" % ( topValue, "\t".join( exampleDict[ topValue]) ) )
            
        outFile.close()
        outFile2.close()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

