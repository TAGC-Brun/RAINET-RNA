import sys
import os
import argparse
from scipy import stats
import cPickle as pickle

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil

#===============================================================================
# Started 20-May-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to read large catRAPID interaction files."
SCRIPT_NAME = "ReadCatrapid.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Parse catRAPID interaction file
#===============================================================================


class ReadCatrapid(object):
    
    STORED_INTERACTIONS_FILENAME = "/storedInteractions_"
    PROTEIN_INTERACTIONS_FILENAME = "/proteinInteractions.tsv"
    
    def __init__(self, catrapid_file, output_folder, interaction_cutoff, interaction_filter_file, batch_size):

        self.catRAPIDFile = catrapid_file
        self.outputFolder = output_folder
        self.interactionCutoff = interaction_cutoff
        self.interactionFilterFile = interaction_filter_file
        self.batchSize = batch_size

    # #
    # Read list of interacting pairs
    # @return set of interacting pairs we want to keep, empty if no file given
    def read_interaction_filter_file(self):
        
        if self.interactionFilterFile != "":
            with open( self.interactionFilterFile, "r") as inFile:        
                wantedPairs = { "_".join( line.strip().split("\t")) for line in inFile}
    
            print "read_interaction_filter_file: read %s unique pairs interacting pairs." % len( wantedPairs)
            
            return wantedPairs
        else:
            return set()
        
    # #
    # Read catrapid file and write processed output to files
    def read_catrapid_file( self, wanted_pairs):

        #=======================================================================
        # Example file
        # sp|Q96DC8|ECHD3_HUMAN ENST00000579524   -12.33  0.10    0.00
        # sp|P10645|CMGA_HUMAN ENST00000516610    10.66   0.32    0.00
        # protein and rna separated by " ", other values separated by "\t"
        # 
        # Protein is always on left side, RNA in the right side.
        # Assumption that there only one interaction between each Protein-RNA pair
        #=======================================================================

        #=======================================================================
        # initialising
        #=======================================================================

        # process interactionCutoff attribute
        if self.interactionCutoff == "OFF":
            self.interactionCutoff = float( "-inf")
        else:
            self.interactionCutoff = float( self.interactionCutoff)

        # check if we need to filter by wanted pairs
        if len( wanted_pairs) > 0:
            filterBool = 1
        else:
            filterBool = 0

        # approach: initialise protein, and sum scores throughout the file, and keep count of protein occurrences, then in the end calculate mean
        proteinInteractions = {} # key -> protein ID, value -> sum of scores
        proteinInteractionsCounter = {} # key -> protein ID, value -> number of times protein appears

        lineCount = 0
        outFileCount = 1

        inFileName = os.path.basename( self.catRAPIDFile)

        # variable which will store text to be written into files        
        interactionText = ""

        #=======================================================================
        # read file
        #=======================================================================
        with open( self.catRAPIDFile, "r") as inFile:
            for line in inFile:

                # every X lines, write to file to liberate memory                
                if lineCount % self.batchSize == 0 and lineCount != 0:
                    Timer.get_instance().step("read_catrapid_file: reading %s lines.." % lineCount)    

                    # print len( proteinInteractions), sys.getsizeof( proteinInteractions) / 1000000.0
                    # print len( interactionText), sys.getsizeof( interactionText) / 1000000.0 

                    # dump dictionaries into files
                    with open( self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME + str( outFileCount) + ".tsv", "w") as outFile:
                        outFile.write( interactionText)
                    interactionText = ""
                    
                    outFileCount += 1

                lineCount += 1 # this has to be before the 'continues'
                    
                spl = line.split(" ")
                
                protID = spl[0].split( "|")[1]
                spl2 = spl[1].split( "\t")
                rnaID = spl2[0]
                score = float( spl2[1])
                
                pair = "_".join( [protID, rnaID])
                                
                # filter by score
                if score < self.interactionCutoff: 
                    continue

                # if filtering by wanted pairs and it is not present
                if filterBool and pair not in wanted_pairs:
                    continue

                # store interaction
                interactionText += "%s\t%s\n" % (pair, score)

                # for calculating average score per protein
                if protID not in proteinInteractions:
                    proteinInteractions[ protID] = 0
                    proteinInteractionsCounter[ protID] = 0
 
                proteinInteractions[ protID] += score
                proteinInteractionsCounter[ protID] += 1


            # write remaining interactions into file
            with open( self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME + str( outFileCount) + ".tsv", "w") as outFile:
                outFile.write( interactionText)

            # write file and dictionary with mean interaction score per protein
            proteinInteractionsMean = {} # key -> protein ID, value -> mean score
            with open( self.outputFolder + ReadCatrapid.PROTEIN_INTERACTIONS_FILENAME, "w") as outFile:
                for prot in proteinInteractions:
                    # sum of scores divided by frequency
                    mean = proteinInteractions[ prot] / float( proteinInteractionsCounter[ prot])
                    proteinInteractionsMean[ prot] = mean
                    outFile.write( "%s\t%s\n" % (prot, mean) )


        print "read_catrapid_file: read %s lines.." % lineCount

        assert len( proteinInteractions) == len( proteinInteractionsCounter) == len( proteinInteractionsMean)

        return proteinInteractionsMean, proteinInteractionsCounter


if __name__ == "__main__":
    
    # Start chrono
    Timer.get_instance().start_chrono()
    
    print "STARTING " + SCRIPT_NAME
    
    #===============================================================================
    # Get input arguments, initialise class
    #===============================================================================
    parser = argparse.ArgumentParser(description= DESC_COMMENT) 

    # positional args
    parser.add_argument('catRAPIDFile', metavar='catRAPIDFile', type=str,
                         help='Output file from catRAPID library all vs all.')
    parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
    parser.add_argument('--interactionCutoff', metavar='interactionCutoff', type=str,
                         default = 15, help='Minimum catRAPID interaction propensity. Set as "OFF" if no filtering wanted.')
    parser.add_argument('--interactionFilterFile', metavar='interactionFilterFile', type=str,
                         default = "", help='TSV file with list of interacting pairs we want to keep, one pair per line. UniprotAC\tEnsemblTxID.')
    parser.add_argument('--batchSize', metavar='batchSize', type=int,
                         default = 1000000, help='How many lines to process before writing to file (to avoid excessive memory consumption).')   

    #gets the arguments
    args = parser.parse_args( ) 

    # init
    run = ReadCatrapid( args.catRAPIDFile, args.outputFolder, args.interactionCutoff, args.interactionFilterFile, args.batchSize)

    # read interaction filter file (if any)
    Timer.get_instance().step( "reading interaction filter file..")    
    wantedPairs = run.read_interaction_filter_file( )

    Timer.get_instance().step( "reading catrapid interaction file..")    
    run.read_catrapid_file( wantedPairs)

    # Stop the chrono      
    Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )
