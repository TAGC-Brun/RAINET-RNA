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
DESC_COMMENT = "Script to read, filter and process large catRAPID interaction files."
SCRIPT_NAME = "ReadCatrapid.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Parse catRAPID interaction file
# 2) apply interaction filters
# 3) write filtered interaction file, and other processed data files
#===============================================================================

#===============================================================================
# Processing notes:
# 1) To reduce memory consumption, the score values are rounded to 1 decimal. 
#    Thus, means are not precise
# 2) Filters are all applied on top of each other, first by score, then RNA, then protein, then interaction-based.
#===============================================================================


class ReadCatrapid(object):
    
    TEMP_STORED_INTERACTIONS_FILENAME = "/temp_storedInteractions_"
    STORED_INTERACTIONS_FILENAME = "/storedInteractions.tsv"
    PROTEIN_INTERACTIONS_FILENAME = "/proteinInteractions.tsv"
    RNA_INTERACTIONS_FILENAME = "/rnaInteractions.tsv"
    ALL_INTERACTIONS_FILTERED_TAG = "NA" # value to give when an RNA or protein has all their interactions filtered with cutoff
    NORMALISED_STORED_INTERACTIONS_FILENAME = "/storedInteractionsNormalised.tsv"
    INTERACTIONS_SCORE_MATRIX = "interaction_score_matrix.tsv"
    MAXIMUM_NUMBER_VIABLE_INTERACTIONS = 10000000 # maximum number of interactions writable for interaction matrix output
    
    def __init__(self, catrapid_file, output_folder, interaction_cutoff, interaction_filter_file, rna_filter_file, protein_filter_file,
                 write_interactions, batch_size, write_normalised_interactions, write_interaction_matrix):

        self.catRAPIDFile = catrapid_file
        self.outputFolder = output_folder
        self.interactionCutoff = interaction_cutoff
        self.interactionFilterFile = interaction_filter_file
        self.rnaFilterFile = rna_filter_file
        self.proteinFilterFile = protein_filter_file
        self.writeInteractions = write_interactions
        self.batchSize = batch_size
        self.writeNormalisedInteractions = write_normalised_interactions
        self.writeInteractionMatrix = write_interaction_matrix

        if (write_normalised_interactions and write_interactions == 0) or (write_interaction_matrix and write_interactions == 0):
            raise RainetException( "ReadCatrapid.__init__ : --write_interactions option must be on for --write_normalised_interactions to run.")

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)
        else:
            print "__init__: Output folder already exists: %s" % self.outputFolder


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
    # Read list of wanted rnas 
    # @return set of rnas we want to keep, empty if no file given
    def read_rna_filter_file(self):
        
        if self.rnaFilterFile != "":
            with open( self.rnaFilterFile, "r") as inFile:        
                wantedList = { line.strip() for line in inFile }
            print "read_interaction_filter_file: read %s unique wanted RNAs." % len( wantedList)           
            return wantedList
        else:
            return set()
        
    # #
    # Read list of wanted proteins
    # @return set of proteins we want to keep, empty if no file given
    def read_protein_filter_file(self):
        
        if self.proteinFilterFile != "":
            with open( self.proteinFilterFile, "r") as inFile:        
                wantedList = { line.strip() for line in inFile }
            print "read_interaction_filter_file: read %s unique wanted Proteins." % len( wantedList)           
            return wantedList
        else:
            return set()
 
    # #
    # Read catrapid file, apply filters, and write processed output to files
    def read_catrapid_file( self, wanted_pairs, wanted_RNAs, wanted_proteins):

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

        # check if we need to filter by wanted pairs, proteins, rnas
        if len( wanted_pairs) > 0: interactionFilterBool = 1
        else: interactionFilterBool = 0

        if len( wanted_RNAs) > 0: rnaFilterBool = 1
        else: rnaFilterBool = 0

        if len( wanted_proteins) > 0: proteinFilterBool = 1
        else: proteinFilterBool = 0

        
        ### Protein containers ####

#         # approach1: initialise protein, and sum scores throughout the file, and keep count of protein occurrences, then in the end calculate mean
#         proteinInteractionsSum = {} # key -> protein ID, value -> sum of scores
        proteinInteractionsCounter = {} # key -> protein ID, value -> number of times protein appears
        proteinInteractionsMean = {}
        # approach2: initialise protein, create a dictionary for each protein which contains the frequencies of each score instead of list of scores, in order to save memory
        proteinScoreFrequencies = {} # key -> protein ID, value -> dict. key -> score, value -> frequency of score
        allProtSet = set()

        ### RNA containers ####
        # approach: initialise RNA, create a dictionary for each RNA which contains the frequencies of each score instead of list of scores, in order to save memory
        rnaScoreFrequencies = {} # key -> RNA ID, value -> dict. key -> score, value -> frequency of score
        allRNASet = set()

        lineCount = 0
        outFileCount = 1

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

                    # print len( proteinInteractionsSum), sys.getsizeof( proteinInteractionsSum) / 1000000.0
                    # print len( interactionText), sys.getsizeof( interactionText) / 1000000.0 

                    # dump dictionaries into files
                    if self.writeInteractions:
                        with open( self.outputFolder + ReadCatrapid.TEMP_STORED_INTERACTIONS_FILENAME + str( outFileCount) + ".tsv", "w") as outFile:
                            outFile.write( interactionText)
                    
                    interactionText = ""
                    
                    outFileCount += 1

                lineCount += 1 # this has to be before the filterings ( 'continue')
                    
                spl = line.split(" ")
                
                protID = spl[0].split( "|")[1]
                spl2 = spl[1].split( "\t")
                rnaID = spl2[0]
                score = float( spl2[1])
                scoreRounded = round( score, 1) 
 
                pair = "_".join( [protID, rnaID])
                           
                allRNASet.add( rnaID)
                allProtSet.add( protID)
               
                #### Apply filterings ####     
                # filter by score
                if score < self.interactionCutoff: 
                    continue

                # if filtering by wanted RNAs and it is not present
                if rnaFilterBool and rnaID not in wanted_RNAs:
                    continue

                # if filtering by wanted Proteins and it is not present
                if proteinFilterBool and protID not in wanted_proteins:
                    continue

                # if filtering by wanted pairs and it is not present
                if interactionFilterBool and pair not in wanted_pairs:
                    continue

                #### Store interaction #### 
                #interactionText += "%s\t%s\n" % (pair, score)
                interactionText+= line

                ## Protein side
#                 # for calculating average score per protein
                if protID not in proteinScoreFrequencies:
                    proteinScoreFrequencies[ protID] = {}

                # producing dictionary with score frequencies for a protein
                if scoreRounded not in proteinScoreFrequencies[ protID]:
                    proteinScoreFrequencies[ protID][ scoreRounded] = 0
                proteinScoreFrequencies[ protID][ scoreRounded] += 1
 
                ## RNA side
                if rnaID not in rnaScoreFrequencies:
                    rnaScoreFrequencies[ rnaID] = {}

                if scoreRounded not in rnaScoreFrequencies[ rnaID]:
                    rnaScoreFrequencies[ rnaID][ scoreRounded] = 0
                rnaScoreFrequencies[ rnaID][ scoreRounded] += 1


            # write remaining interactions into file
            if self.writeInteractions:
                with open( self.outputFolder + ReadCatrapid.TEMP_STORED_INTERACTIONS_FILENAME + str( outFileCount) + ".tsv", "w") as outFile:
                    outFile.write( interactionText)

        print "read_catrapid_file: read %s lines.." % lineCount

        # join stored interaction files
        if self.writeInteractions:
            # cat files
            cmd = "cat %s* > %s" % ( self.outputFolder + ReadCatrapid.TEMP_STORED_INTERACTIONS_FILENAME, self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME )
            os.system( cmd)
            
            # remove temp files
            cmd = "rm %s*" % ( self.outputFolder + ReadCatrapid.TEMP_STORED_INTERACTIONS_FILENAME)
            os.system( cmd)


        #=======================================================================
        # Write output file
        #=======================================================================

        ### RNA file ###
 
        with open( self.outputFolder + ReadCatrapid.RNA_INTERACTIONS_FILENAME, "w") as outFile:
            # change header here
            outFile.write("ensembl_id\tmean_score\tmedian_score\tmin_score\tmax_score\tstd_score\tcount\n")

            for rna in allRNASet:
                if rna in rnaScoreFrequencies:
                    # recreate all original values for a rna
                    listOfScores = [ scoreVal for scoreVal in rnaScoreFrequencies[ rna] for i in range( rnaScoreFrequencies[ rna][ scoreVal])]
     
                    mean = np.mean( listOfScores)
                    median = np.median( listOfScores)   
                    minimum = np.min( listOfScores)
                    maximum = np.max( listOfScores)
                    std = np.std( listOfScores)    
                    # number of Proteins/interactions above filter
                    count = len(listOfScores)
                    
                    outFile.write( "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n" % (rna, mean, median, minimum, maximum, std, count) )

                else:                   
                    outFile.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( rna, 
                                                              ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG, ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG,
                                                              ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG, ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG,
                                                              ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG, ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG ) )


        ### Protein file ###

        with open( self.outputFolder + ReadCatrapid.PROTEIN_INTERACTIONS_FILENAME, "w") as outFile:
            # change header here
            outFile.write("uniprotac\tmean_score\tmedian_score\tmin_score\tmax_score\tstd_score\tcount\n")

            # calculate protein score metrics
            for prot in allProtSet:
                if prot in proteinScoreFrequencies:

                    # recreate all original values for a protein
                    listOfScores = [ scoreVal for scoreVal in proteinScoreFrequencies[ prot] for i in range( proteinScoreFrequencies[ prot][ scoreVal])]
    
                    mean = np.mean( listOfScores)
                    median = np.median( listOfScores)
                    minimum = np.min( listOfScores)
                    maximum = np.max( listOfScores)
                    std = np.std( listOfScores)
                        
                    # number of RNAs above filter
                    count = len(listOfScores)
    
                    proteinInteractionsMean[ prot] = mean
                    proteinInteractionsCounter[ prot] = count
                    
                    outFile.write( "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n" % (prot, mean, median, minimum, maximum, std, count) )
    
                else:                   
                    outFile.write( "%s\t%s\t%s\t\%st%s\t%s\t%s\n" % ( prot, 
                                                              ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG, ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG,
                                                              ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG, ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG,
                                                              ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG, ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG ) )

        return proteinInteractionsMean, proteinInteractionsCounter


    # #
    # Function to write extra output file with interactions normalised by the max for each RNA, 
    # (using unity-based normalisation, aka min-max normalisation)
    # This function runs after writing non-normalised interactions file and rna interactions file.
    def write_normalised_interactions( self):

        if not os.path.exists( self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME):
            raise RainetException( "ReadCatrapid.write_normalised_interactions : output interactions file not found. %s" % self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME)
        if not os.path.exists( self.outputFolder + ReadCatrapid.RNA_INTERACTIONS_FILENAME):
            raise RainetException( "ReadCatrapid.write_normalised_interactions : output rna interactions file not found. %s" % self.outputFolder + ReadCatrapid.RNA_INTERACTIONS_FILENAME)


        #===============================================================================
        # Get maximum and minimum scores for each transcript
        #===============================================================================

        # e.g. format: ensembl_id      mean_score      median_score    min_score       max_score       std_score       count
        #        ENST00000388090 -15.11  -16.10  -45.90  26.90   10.46   1978

        # column indexes
        txIDCol = 0
        minCol = 3 
        maxCol = 4

        rnaMax = {} # key -> transcriptID, val -> max score among all interactions
        rnaMin = {} # key -> transcriptID, val -> min score among all interactions
        with open( self.outputFolder + ReadCatrapid.RNA_INTERACTIONS_FILENAME, "r") as inFile:
            inFile.readline() # skip header
            
            for line in inFile:
                line = line.strip()
                spl = line.split( "\t")
                
                transcriptID = spl[ txIDCol]
                minimum = spl[ minCol]
                maximum = spl[ maxCol]

                # if this RNA was filtered out, it will also not feature in the interactions file.
                if minimum == ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG or maximum == ReadCatrapid.ALL_INTERACTIONS_FILTERED_TAG:
                    continue

                minimum = float( minimum)
                maximum = float( maximum)

                if transcriptID not in rnaMax:
                    rnaMax[ transcriptID] = maximum
                else:
                    raise RainetException( "ReadCatrapid.write_normalised_interactions : duplicate transcript ID. %s" % transcriptID)

                if transcriptID not in rnaMin:
                    rnaMin[ transcriptID] = minimum
                else:
                    raise RainetException( "ReadCatrapid.write_normalised_interactions : duplicate transcript ID. %s" % transcriptID)


        #===============================================================================
        # Apply normalisation
        #===============================================================================

        idColumn = 0        
        scoreColumn = 1
        
        # e.g. format: sp|Q7Z419|R144B_HUMAN ENST00000542804    20.56    0.54    0.00

        outFile = open( self.outputFolder + ReadCatrapid.NORMALISED_STORED_INTERACTIONS_FILENAME, "w")
        
        with open( self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME, "r") as inFile:
            
            for line in inFile:
                line = line.strip()
                spl = line.split( "\t")
                
                score = round( float( spl[ scoreColumn]), 1)

                transcriptID = spl[ idColumn].split( " ")[1]
                                
                # min-max normalisation

                minimum = rnaMin[ transcriptID]
                maximum = rnaMax[ transcriptID]
                
                assert score >= minimum
                assert score <= maximum

                normalisedScore = self._min_max_normalisation( score, minimum, maximum)
                
                # rewrite output file
                text = "%s\t%.2f\t%s\n" % ( "\t".join( spl[ :scoreColumn]), normalisedScore, "\t".join( spl[ scoreColumn+1:]))

                outFile.write( text)

        outFile.close()
        

    # function to calculate min-max normalisation (unity-based normalisation)
    def _min_max_normalisation(self, x, minimum, maximum):

        normVal = float( x - minimum) / float( maximum - minimum)
       
        assert normVal >= 0        
        assert normVal <= 1

        return normVal


    # #
    # Function to write matrix output file for interactions after filtering.
    # This function runs after writing interactions file.
    def write_matrix_output( self):

        #=================================================================== 
        # Initialising
        #=================================================================== 

        if not os.path.exists( self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME):
            raise RainetException( "ReadCatrapid.write_matrix_output : output interactions file not found. %s" % self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME)

        ## Check amount of interactions to not overload the system
        
        cmd = "wc -l %s" % ( self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME)

        stout = SubprocessUtil.run_command( cmd, verbose = 0, return_stdout = 1)
        try:
            nlines = int( stout.split( " ")[0])
        except:
            raise RainetException( "Could not calculate number of lines in filtered stored interactions." % cmd)

        if nlines > ReadCatrapid.MAXIMUM_NUMBER_VIABLE_INTERACTIONS:
            raise RainetException( "ReadCatrapid.write_matrix_output : number of interactions to write matrix is too large to be computable: %s interactions" % nlines)
            
        print "Writing matrix with %s interactions" % nlines

        #=================================================================== 
        # Read filtered interactions file and store interactions into memory
        #=================================================================== 

        # create data structures with all proteins, RNAs and scores of pairs 
        setInteractingRNAs = set()
        setInteractingProts = set()
        dictPairs = {}

        idColumn = 0        
        scoreColumn = 1
         
        # e.g. format: sp|Q7Z419|R144B_HUMAN ENST00000542804    20.56    0.54    0.00
  
        with open( self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME, "r") as inFile:
             
            for line in inFile:
                line = line.strip()
                spl = line.split( "\t")
                 
                score = float( spl[ scoreColumn])
 
                spl2 = spl[ idColumn].split( " ")
 
                txID = spl2[ 1]
                protID = spl2[ 0]
                            
                pair = txID + "|" + protID

                if pair not in dictPairs:
                    dictPairs[pair] = score
                else:
                    raise RainetException("ReadCatrapid.write_matrix_output: duplicate interaction", pair)
     
                setInteractingRNAs.add( txID)
                setInteractingProts.add( protID)
 
        #=================================================================== 
        # Write file with interaction scores for each protein-RNA pair, matrix format
        #=================================================================== 
        # E.g.
        # RNAs Prot1 Prot2
        # RNA1 10.4    0.3
        # RNA2 32.6    -34.5
 
        outHandler = open( self.outputFolder + ReadCatrapid.INTERACTIONS_SCORE_MATRIX, "w")
 
        # use sorting to keep headers in place
        sortedSetInteractingProts = sorted( setInteractingProts)
        sortedSetInteractingRNAs = sorted( setInteractingRNAs)
 
        # write header with protein IDs
        outHandler.write( "RNAs")
        for prot in sortedSetInteractingProts:
            outHandler.write( "\t%s" % prot )
        outHandler.write( "\n")
             
        # write bulk of file, one row per rna, one column per protein
        for rna in sortedSetInteractingRNAs:
            text = rna
            for prot in sortedSetInteractingProts:
                tag = rna + "|" + prot
                if tag in dictPairs:
                    score = dictPairs[tag]
                else:
                    score = "NA"
                text+= "\t%s" % score 
            text+= "\n"
            outHandler.write( text)
 
        outHandler.close()

        
    # run functions in proper order
    def run(self):
        
        Timer.get_instance().step( "reading filter files..")    
        wantedPairs = self.read_interaction_filter_file( )
        wantedRNAs = self.read_rna_filter_file( )
        wantedProteins = self.read_protein_filter_file( )
    
        Timer.get_instance().step( "reading catrapid interaction file..")    
        self.read_catrapid_file( wantedPairs, wantedRNAs, wantedProteins)

        if self.writeNormalisedInteractions:
            Timer.get_instance().step( "writing normalised interaction file..")    
            self.write_normalised_interactions( )

        if self.writeInteractionMatrix:
            Timer.get_instance().step( "writing interaction matrix file..")    
            self.write_matrix_output()
            


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
                         default = "OFF", help='Minimum catRAPID interaction propensity. Set as "OFF" if no filtering wanted.')
    parser.add_argument('--interactionFilterFile', metavar='interactionFilterFile', type=str,
                         default = "", help='TSV file with list of interacting pairs we want to keep, one pair per line. UniprotAC\tEnsemblTxID. No header.')
    parser.add_argument('--rnaFilterFile', metavar='rnaFilterFile', type=str,
                         default = "", help='File with list of RNAs we want to keep, one per line. No header.')
    parser.add_argument('--proteinFilterFile', metavar='proteinFilterFile', type=str,
                         default = "", help='File with list of Proteins we want to keep, one per line. No header.')
    parser.add_argument('--writeInteractions', metavar='writeInteractions', type=int,
                         default = 1, help='Whether to write interaction file after the filters.')
    parser.add_argument('--batchSize', metavar='batchSize', type=int,
                         default = 1000000, help='How many lines to process before writing to file (to avoid excessive memory consumption).')   
    parser.add_argument('--writeNormalisedInteractions', metavar='writeNormalisedInteractions', type=int,
                         default = 0, help='Whether to write interaction file after the filters, normalised by max (unity-based normalisation) score for each RNA. --writeInteractions argument must also be 1.')   
    parser.add_argument('--writeInteractionMatrix', metavar='writeInteractionMatrix', type=int,
                         default = 0, help='Whether to write interaction matrix file after the filters. --writeInteractions argument must also be 1.')   

    #gets the arguments
    args = parser.parse_args( ) 

    # init
    readCatrapid = ReadCatrapid( args.catRAPIDFile, args.outputFolder, args.interactionCutoff, args.interactionFilterFile, 
                                 args.rnaFilterFile, args.proteinFilterFile, args.writeInteractions, args.batchSize, 
                                 args.writeNormalisedInteractions, args.writeInteractionMatrix)

    readCatrapid.run()

    # Stop the chrono      
    Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

