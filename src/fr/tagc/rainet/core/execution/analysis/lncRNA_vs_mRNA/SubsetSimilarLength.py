
import sys
import os
import argparse
import numpy as np
from scipy import stats
import random

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

#===============================================================================
# Started 14-June-2016
# Diogo Ribeiro
# Script to create a subset of mRNAs that have similar length distribution as lncRNAs.
# Objective is to control for length of transcript when comparing scores
DESC_COMMENT = "Script to create a subset of mRNAs that have similar length distribution as lncRNAs."
SCRIPT_NAME = "SubsetSimilarLength.py"
#===============================================================================


class SubsetSimilarLength(object):

    MINIMUM_LENGTH = 50
    
    def __init__(self, lnc_rnalengths, m_rnalengths, output_folder, bin_size, length_max):

        self.lncRNAlengths = lnc_rnalengths
        self.mRNAlengths = m_rnalengths
        self.outputFolder = output_folder
        self.binSize = bin_size
        self.lengthMax = length_max

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)

    # #
    # Function to put transcripts into length bins
    def read_fastalength(self, in_file_path):

        print "Input file:", in_file_path

        ### create bins of lengths
        bins = {} # key -> binName, val -> list of transcripts in bin
        for i in range(SubsetSimilarLength.MINIMUM_LENGTH, self.lengthMax, self.binSize):
            binName = str(i) + "-" + str( min( i + self.binSize, self.lengthMax))
            
            bins[ binName] = set()

        print "Number bins created:", len( bins)
        
        ### read lengths file and put transcripts into bins
        # Format e.g. (space separated)
        # length transcript
        # 915 ENST00000510508
        notBinned = 0
        countLines = 0
        lengthDict = {} # key -> transcript, value -> length
        with open( in_file_path, "r") as inFile:
            inFile.readline() #skip header
            for line in inFile:
                countLines += 1
                line = line.strip()                
                spl = line.split(" ")
                length = int( spl[0])
                txID = spl[1]
                
                lengthDict[ txID] = length
                
                boo = 0

                for binName in bins:
                    binSpl = binName.split("-")
                    binStart = int( binSpl[0])
                    binEnd = int( binSpl[1])

                    if length >= binStart and length < binEnd:
                        boo = 1
                        bins[ binName].add( txID)

                if boo == 0:
                    notBinned += 1
#                    print "Could not find bin for transcript:", line
                    
        countTx = 0
        for binName in bins:
#            print binName, bins[binName]
            countTx += len(bins[ binName])

        print "Transcripts binned:", countTx
        print "Transcripts not binned:", notBinned

        assert countLines == countTx + notBinned

        return bins, lengthDict


    # #
    # Given two distributions in bins, make two lists of items that will have similar distribution
    # Priority is given for bins2 number of items to match bins1, but bins1 will also match bins2 in case bins2 has less items.
    def match_distributions(self, bins1, bins2, lengths1, lengths2):

        newList1 = set()
        newList2 = set()
        
        for binName in bins1:
                        
            # transcript lists
            tBins1 = bins1[ binName]
            tBins2 = bins2[ binName]

            # number of transcripts in bin
            nBins1 = len( bins1[ binName])
            nBins2 = len( bins2[ binName])

            if nBins1 > nBins2:
                # downsample nBins1 to match nBins2
                subset = random.sample( tBins1, nBins2)
                newList1.update( set( subset))
                newList2.update( tBins2)             
            elif nBins2 > nBins1:
                # downsample nBins2 to match nBins1
                subset = random.sample( tBins2, nBins1)
                newList1.update( tBins1)
                newList2.update( set( subset))
            else:
                # equality
                newList1.update( tBins1)
                newList2.update( tBins2)


        # confirm that two new lists have same number of transcripts
        assert len(newList1) == len(newList2)

        print "Number of transcripts on each list:", len(newList1)

        ## write to output file
        
        with open( self.outputFolder + "/list_lncRNAs.txt","w") as outFile:
            outFile.write("ensembl_id\tlength\n")
            for item in newList1:
                outFile.write( "%s\t%s\n" % (item, lengths1[item]))
        with open( self.outputFolder + "/list_mRNAs.txt","w") as outFile:
            outFile.write("ensembl_id\tlength\n")
            for item in newList2:
                outFile.write( "%s\t%s\n" % (item, lengths2[item]))
        

if __name__ == "__main__":
    
    # Start chrono
    Timer.get_instance().start_chrono()
    
    print "STARTING " + SCRIPT_NAME
    
    #===============================================================================
    # Get input arguments, initialise class
    #===============================================================================
    parser = argparse.ArgumentParser(description= DESC_COMMENT) 

    # positional args
    parser.add_argument('lncRNAlengths', metavar='lncRNAlengths', type=str,
                         help='Output from fastalength')
    parser.add_argument('mRNAlengths', metavar='mRNAlengths', type=str,
                         help='Output from fastalength')
    parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
    parser.add_argument('--binSize', metavar='binSize', default = 10, type=int, help='Length bin size.')
    parser.add_argument('--lengthMax', metavar='lengthMax', default = 1200, type=int, help='The maximum length size for both libraries.')


    #gets the arguments
    args = parser.parse_args( ) 

    # init
    run = SubsetSimilarLength( args.lncRNAlengths, args.mRNAlengths, args.outputFolder, args.binSize, args.lengthMax)
    
    lncRNABins, lncRNALengths = run.read_fastalength( args.lncRNAlengths)
    mRNABins, mRNALengths = run.read_fastalength( args.mRNAlengths)
    
    run.match_distributions( lncRNABins, mRNABins, lncRNALengths, mRNALengths)
    
    # Stop the chrono      
    Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )


