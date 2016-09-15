import sys
import os
import argparse
import re
import numpy as np

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

#===============================================================================
# Started 15-Sep-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to combine results of catRAPID all vs all on different fragments of same transcript. Merge of scores of fragments from same RNA can be done by mean or max score."
SCRIPT_NAME = "combine_catrapid_fragment_results.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read (special) interactions file, store all results between same transcript (but different fragments) and same protein.
# 2) Use either mean or max score to transform results from several fragments to a single value. Rewrite interactions file.
#===============================================================================

#===============================================================================
# Processing notes:
# 1)
# 2)
#===============================================================================

# #
# Reads interaction file, combines fragment information and writes output interaction file
def read_interactions( interactions_file, output_file, use_mean):


    ########################################################
    # read interactions file, store information
    ########################################################

    # E.g. of "fragments" interaction file
    # sp|Q5VZY2|PLPP4_HUMAN 1ENST00000423456_MEG3_002.fa_1650-1770    20.80    0.54    0.00
    # sp|Q5VZY2|PLPP4_HUMAN 1ENST00000423456_MEG3_002.fa_1653-1765    19.47    0.52    0.00

    pairs = {} # key -> tag of protein-RNA, value -> list of interaction scores for the several fragments

    tagSeparator = "/"

    with open( interactions_file, "r") as inFile:
                
        for line in inFile:
            line = line.strip()

            spl = line.split("\t")

            idLine = spl[0]
            score = float( spl[1])
            
            # retrieve protein ID and base transcript ID (ENST*)
            protID, rnaID = idLine.split(" ")
            
            # retrieve only the actually ENST ID
            query = re.search(".*(ENST[0-9]+)_?" , rnaID)
            
            if query != None:
                enstID = query.group(1)
            else:
                raise RainetException( "No ENST ID found in line: " + line)

            tag = protID + tagSeparator + enstID
            
            if tag not in pairs:
                pairs[ tag] = []
            
            pairs[ tag].append( score)

    # quick validation
    # 15974 * 3 = 47922
    # print len(pairs) --> 47922
    # there was 3 transcripts as input, and 60 fragments, when combined they give 47922 total interactions

    ########################################################
    # merge values and write output file    
    ########################################################
    
    outFile = open( output_file, "w")
    
    for tag in pairs:
        
        if use_mean:
            newScore = "%.2f" % np.mean( pairs[tag])
        else:
            newScore = "%.2f" % np.max( pairs[tag])

        spl = tag.split( tagSeparator)
        originalTag = spl[0] + " " + spl[1]
        rest = [originalTag, newScore, "NA", "NA"]
        text = "\t".join( rest) + "\n"

        outFile.write( text)

    outFile.close()

    # quick validation
    # sp|Q8WY07|CTR3_HUMAN/ENST00000423456 19.45
    # sp|Q8WY07|CTR3_HUMAN/ENST00000423456 0.781764705882
    # got same results with manual calculation


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
        parser.add_argument('interactionsFile', metavar='interactionsFile', type=str,
                             help='CatRAPID (omics) all vs all interactions file for fragmented transcripts. ENST* ID must exist for the RNA. Example format: sp|Q5VZY2|PLPP4_HUMAN 1ENST00000423456_MEG3_002.fa_1650-1770    20.80    0.54    0.00')
        parser.add_argument('outputFile', metavar='outputFile', type=str,
                             help='Output interactions files with merged interaction scores between fragments. Other score metrics are not passed to output.')
        parser.add_argument('--useMean', metavar='useMean', type=int, default = 0,
                             help='If turned on, use mean instead of max for merging results from several fragments of same transcript (Default = 0).')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        # run function to..
        Timer.get_instance().step( "Read catRAPID file file..")            
        read_interactions( args.interactionsFile, args.outputFile, args.useMean)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

