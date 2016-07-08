import sys
import os
import argparse
import numpy as np
from scipy import stats
import cPickle as pickle

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil

#===============================================================================
# Started 20-May-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to read Starbase file and create file with line per protein, as in ReadCatrapid.py."
SCRIPT_NAME = "read_starbase.py"
#===============================================================================

#===============================================================================
# General plan:
#===============================================================================

#===============================================================================
# Processing notes:
#===============================================================================

# #
# Read starbase file, write file with measures per protein
def read_starbase_file( star_base_file, output_folder, minimum_bio_complex, minimum_clip_read_number):

    #=======================================================================
    # Example file
    # pairID  clipReads       bioComplex
    # ENST00000428677|Q8TBR3  18      2    # 
    # RNA is always on left side, Protein in the right side.
    #=======================================================================

    #=======================================================================
    # initialising
    #=======================================================================
    
    ### Protein containers ####

    # approach2: initialise protein, create a dictionary for each protein which contains the frequencies of each score instead of list of scores, in order to save memory
    proteinScoreFrequencies = {} # key -> protein ID, value -> dict. key -> score, value -> frequency of score
    allProtSet = set()

    lineCount = 0

    #=======================================================================
    # read file
    #=======================================================================
    with open( star_base_file, "r") as inFile:
        inFile.readline()
        for line in inFile:
                        
            line = line.strip()

            lineCount += 1
                
            spl = line.split("\t")
            
            rnaID, protID = spl[0].split( "|")
            clipReads = float( spl[1])
            bioComplex = float( spl[2])

            # apply filterings            
            if clipReads < minimum_clip_read_number:
                continue
            if bioComplex < minimum_bio_complex:
                continue
                       
            allProtSet.add( protID)
           
            ## Protein side
#                 # for calculating average score per protein
            if protID not in proteinScoreFrequencies:
                proteinScoreFrequencies[ protID] = {}

            # producing dictionary with score frequencies for a protein
            if clipReads not in proteinScoreFrequencies[ protID]:
                proteinScoreFrequencies[ protID][ clipReads] = 0
            proteinScoreFrequencies[ protID][ clipReads] += 1


    print "read_starbase_file: read %s lines.." % lineCount

    #=======================================================================
    # Write output file
    #=======================================================================


    ### Protein file ###

    with open( output_folder + "/protein_interactions.tsv", "w") as outFile:
        # change header here
        outFile.write("uniprotac\tmean_reads\tmedian_reads\tstd_score\tcount_interactions\n")

        # calculate protein score metrics
        for prot in allProtSet:
            if prot not in proteinScoreFrequencies:
                print "PROBLEM"

            # recreate all original values for a protein
            listOfScores = [ scoreVal for scoreVal in proteinScoreFrequencies[ prot] for i in range( proteinScoreFrequencies[ prot][ scoreVal])]

            mean = np.mean( listOfScores)
            median = np.median( listOfScores)   
            std = np.std( listOfScores)
                
            # number of RNAs
            count = len(listOfScores)

            outFile.write( "%s\t%.2f\t%.2f\t%.2f\t%s\n" % (prot, mean, median, std, count) )




if __name__ == "__main__":
    
    # Start chrono
    Timer.get_instance().start_chrono()
    
    print "STARTING " + SCRIPT_NAME
    
    #===============================================================================
    # Get input arguments, initialise class
    #===============================================================================
    parser = argparse.ArgumentParser(description= DESC_COMMENT) 

    # positional args
    parser.add_argument('starBaseFile', metavar='starBaseFile', type=str,
                         help='starbase_interactions.tsv from StarBasePredictionValidation.py.')
    parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
    parser.add_argument('--minimumBioComplex', metavar='minimumBioComplex', default = 0, type=int, help='Minimum value of BioComplex to keep starBase interaction.')
    parser.add_argument('--minimumClipReadNumber', metavar='minimumClipReadNumber', default = 0, type=int, help='Minimum value of clipReadNum to keep starBase interaction.')

    #gets the arguments
    args = parser.parse_args( ) 

    # init
    read_starbase_file( args.starBaseFile, args.outputFolder, args.minimumBioComplex, args.minimumClipReadNumber)

    # Stop the chrono      
    Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )
