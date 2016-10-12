
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
DESC_COMMENT = "Script to correlate catRAPID scores between two files sharing common interaction pairs"
SCRIPT_NAME = "score_correlation.py"
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
# Read catRAPID file.
# Updated for new catRAPID format. No need for cross references.
def read_catrapid_file_new( input_file):

    # E.g.: sp|Q6P6C2|ALKB5_HUMAN ENST00000559683   47.85   0.93    0.23

    interactingPairs = {} # key -> pair of transcriptID and proteinID, val -> score
    
    proteinSet = set()

    countLines = 0
    
    with open( input_file, "r") as f:
        for line in f:
            spl = line.split(" ")

            countLines+= 1 

            if countLines % 10000000 == 0:
                print "Processed %s interactions" % countLines
                           
            
            proteinID = spl[0].split( "|")[1]
            spl2 = spl[1].split( "\t")
            transcriptID = spl2[0]
            intScore = float( spl2[1])
            
            pair = transcriptID + "|" + proteinID

            proteinSet.add(proteinID)

            # add pair to interacting pairs and keep the maximum interaction score
            if pair not in interactingPairs:
                interactingPairs[ pair] = intScore
            else:
                raise RainetException( "Repeated protein-RNA pair: " + line)


    print "read_catrapid_file_new: Number of proteins: ", len( proteinSet)
    print "read_catrapid_file_new: Number of protein-RNA pairs in catRAPID: ", len( interactingPairs)

    return interactingPairs


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
        parser.add_argument('catrapidFile1', metavar='catrapidFile1', type=str,
                             help='Description.')
        parser.add_argument('catrapidFile2', metavar='catrapidFile2', type=str,
                             help='Description.')
        parser.add_argument('outputFile', metavar='outputFile', type=str,
                             help='Description.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        Timer.get_instance().step( "Read file 1..")            
        pairs1 = read_catrapid_file_new( args.catrapidFile1)
        Timer.get_instance().step( "Read file 2..")            
        pairs2 = read_catrapid_file_new( args.catrapidFile2)

        commonPairs = set( pairs1.keys()).intersection( pairs2.keys())

        print "Number of common pairs: %s" % len( commonPairs)

        outputFile = open( args.outputFile, "w")
        outputFile.write("Pair\tScoreFile1\tScoreFile2\n")

        for pair in commonPairs:
            outputFile.write( "%s\t%s\t%s\n" % (pair, pairs1[ pair], pairs2[ pair]) )

        outputFile.close()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

