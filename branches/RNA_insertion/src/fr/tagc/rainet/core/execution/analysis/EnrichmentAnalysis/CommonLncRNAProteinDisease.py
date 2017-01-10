
import sys
import os
import argparse

# import numpy as np
# import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

# from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
# from fr.tagc.rainet.core.util.data.DataManager import DataManager


#===============================================================================
# Started 05-January-2017
# Diogo Ribeiro
DESC_COMMENT = "Script to match disease from lncRNA to diseases associated to proteins enriched as their respective lncRNA targets."
SCRIPT_NAME = "CommonLncRNAProteinDisease.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read lncRNA disease association file
# 2) Read protein disease association file
# 3) Read lncRNA-protein association file (enrichment)
# 4) Return file with their correspondances, doing a simple word matching
#===============================================================================

#===============================================================================
# Processing notes:
# 1) For the word matching, strings are lowercased and each transcript disease word tries to find a match in the protein disease string using "contains" approach
#===============================================================================

class CommonLncRNAProteinDisease(object):

    # Constants

    OUTPUT_FILE = "/lncRNA_protein_disease_descriptions.txt"
    OUTPUT_FILE_WORD_MATCH = "/lncRNA_protein_disease_descriptions_word_match.txt"

    def __init__(self, lncRNADiseaseFile, proteinDiseaseFile, lncRNAProteinFile, outputFolder, minWordSize, blackListedWords):

        self.lncRNADiseaseFile = lncRNADiseaseFile
        self.proteinDiseaseFile = proteinDiseaseFile
        self.lncRNAProteinFile = lncRNAProteinFile
        self.outputFolder = outputFolder
        self.minWordSize = minWordSize
        self.blackListedWords = blackListedWords

        # stores words to be black listed when performing word matching
        self.blackListedWordsSet = set()

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)
     
            
    # #
    # Read list of words to be black listed when performing word matching
    def read_black_listed_words_file(self):
        
        with open( self.blackListedWords, "r") as inFile:
            for line in inFile:
                self.blackListedWordsSet.add( line.strip())
           
    # #
    # Read lncRNA-disease associations
    def read_lncrna_disease_file(self):

        #===============================================================================
        # Read lncRNA disease file
        #===============================================================================
        # Example format
        # ENST00000475441 endometrial carcinoma
        # ENST00000579368 gastric cancer

        # Note: there is one entry for each transcriptID-disease pair, a transcriptID may have several disease associations and vice-versa
        # Note: file produced manually

        transcriptDiseaseDict = {} # key -> transcriptID, value -> disease name/description
        transcriptDiseases = set()
        nAssociations = 0

        with open( self.lncRNADiseaseFile, "r") as inFile:
            
            # no header
            
            for line in inFile:
                spl = line.strip().split( "\t")
                
                transcriptID = spl[0]
                disease = spl[1]
                
                if transcriptID not in transcriptDiseaseDict:
                    transcriptDiseaseDict[ transcriptID] = set()
                    
                transcriptDiseaseDict[ transcriptID].add( disease)
                transcriptDiseases.add( disease)
                nAssociations += 1

        print "read_lncrna_disease_file: %s transcriptID with disease" % len( transcriptDiseaseDict)
        print "read_lncrna_disease_file: %s transcript diseases" % len( transcriptDiseases)
        print "read_lncrna_disease_file: %s transcript-disease associations" % nAssociations

        self.transcriptDiseaseDict = transcriptDiseaseDict

    # #
    # Read protein-disease associations
    def read_protein_disease_file(self):

        #===============================================================================
        # Read protein disease file
        #===============================================================================
        # Example format
        # O00623  266510  PEROXISOME BIOGENESIS DISORDER 3B; PBD3B |  | 
        # O00623  614859  PEROXISOME BIOGENESIS DISORDER 3A (ZELLWEGER); PBD3A |  | PEROXISOME BIOGENESIS DISORDER, COMPLEMENTATION GROUP 3, INCLUDED; CG3, INCLUDED

        # Note: there is one entry for each proteinID-disease pair, a proteinID may have several disease associations and vice-versa
        # Note: file produced with OMIMProteinDisease.py
    

        proteinDiseaseDict = {} # key -> proteinID, value -> disease name/description
        protenDiseases = set()
        nAssociations = 0

        with open( self.proteinDiseaseFile, "r") as inFile:
            
            # no header
            
            for line in inFile:
                spl = line.strip().split( "\t")
                
                proteinID = spl[0]
                disease = spl[2]
                
                if proteinID not in proteinDiseaseDict:
                    proteinDiseaseDict[ proteinID] = set()
                    
                proteinDiseaseDict[ proteinID].add( disease)
                protenDiseases.add( disease)
                nAssociations += 1

        print "read_protein_disease_file: %s proteinID with disease" % len( proteinDiseaseDict)
        print "read_protein_disease_file: %s protein diseases" % len( protenDiseases)
        print "read_protein_disease_file: %s protein-disease associations" % nAssociations

        self.proteinDiseaseDict = proteinDiseaseDict


    # #
    # Read transcript-protein enrichment correspondence and write output file with disease associations
    def read_lncrna_protein_file(self):

        print "read_lncrna_protein_file: avoiding word match of words below size %s." % self.minWordSize
        
        #===============================================================================
        # Read protein disease file
        #===============================================================================
        # Example format
        # ENST00000424191 P04899
        # ENST00000424191 P09471
        
        # File with original disease terms
        outFile = open( self.outputFolder + CommonLncRNAProteinDisease.OUTPUT_FILE, "w")
        outFile.write("transcriptID\tproteinID\ttranscriptDisease\tproteinDisease\n")

        # File with matched disease terms only
        outFile2 = open( self.outputFolder + CommonLncRNAProteinDisease.OUTPUT_FILE_WORD_MATCH, "w")        
        outFile2.write( "transcriptID\tproteinID\twordsMatched\ttranscriptDisease\tproteinDisease\n")

        # line counters after writing header
        nLines = 1
        nLinesMatch = 1
        
        with open( self.lncRNAProteinFile, "r") as inFile:
            
            # no header
            
            for line in inFile:
                spl = line.strip().split( "\t")
        
                transcriptID = spl[0]
                proteinID = spl[1]
                
                for transcriptDisease in self.transcriptDiseaseDict[ transcriptID]:

                    ## Word matching
                    # put all text in lowercase for word matching
                    txDiseaseLower = transcriptDisease.lower()
                    # Split per word
                    wordSoup = txDiseaseLower.split( " ")
                    # Filter out small words
                    wordSoup = [word for word in wordSoup if len( word) > self.minWordSize]
                                        
                    if proteinID in self.proteinDiseaseDict:
                        for proteinDisease in self.proteinDiseaseDict[ proteinID]:

                            # put text in lowercase to match transcript diseases
                            protDiseaseLower = proteinDisease.lower()

                            # loop each transcript disease word and try to find a match                            
                            wordBoo = 0
                            matchWord = []
                            for word in wordSoup:
                                if word in protDiseaseLower:
                                    wordBoo = 1
                                    # if several words are matched, they are separated by commas (as we dont look at word order) 
                                    matchWord.append( word )

                            if wordBoo:
                                matchWord = ";".join( matchWord).strip()
                                
                                if matchWord not in self.blackListedWordsSet:
                                    outFile2.write( "%s\t%s\t'%s'\t%s\t%s\n" % ( transcriptID, proteinID, matchWord, transcriptDisease, proteinDisease) )
                                    nLinesMatch += 1
                            
                            outFile.write( "%s\t%s\t%s\t%s\n" % ( transcriptID, proteinID, transcriptDisease, proteinDisease) )
                            
                            nLines += 1

        outFile.close()
        outFile2.close()
        
        print "read_lncrna_protein_file: wrote %s association lines." % nLines
        print "read_lncrna_protein_file: wrote %s association lines with a word match." % nLinesMatch


if __name__ == "__main__":

    try:
    
        # Start chrono
        Timer.get_instance().start_chrono()
        
        print "STARTING " + SCRIPT_NAME
        
        #===============================================================================
        # Get input arguments, initialise class
        #===============================================================================
        parser = argparse.ArgumentParser(description= DESC_COMMENT) 
    
        # positional args
        parser.add_argument('lncRNADiseaseFile', metavar='lncRNADiseaseFile', type=str,
                             help='TSV file with lncRNA-disease associations. One line per pair. E.g. ENST00000428939 colorectal cancer')
        parser.add_argument('proteinDiseaseFile', metavar='proteinDiseaseFile', type=str,
                             help='TSV file with protein-disease associations (e.g. from OMIM). One line per pair. E.g. O00623  266510  PEROXISOME BIOGENESIS DISORDER 3B; PBD3B |  | ')
        parser.add_argument('lncRNAProteinFile', metavar='lncRNAProteinFile', type=str,
                             help='TSV file with lncRNA-protein associations, based on lncRNA enrichment to a complex. One line per pair. E.g. ENST00000424191 P63218')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
        parser.add_argument('--minWordSize', metavar='minWordSize', type=str, default = 2,
                             help='For matching transcript disease to protein disease, try to match any words above size X.')
        parser.add_argument('--blackListedWords', metavar='blackListedWords', type=str, default = "",
                             help='Optional file with list of words (one per line) that should not be matched by themselves, unless in conjunction with other words.')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        instance = CommonLncRNAProteinDisease( args.lncRNADiseaseFile, args.proteinDiseaseFile, args.lncRNAProteinFile, args.outputFolder, args.minWordSize, args.blackListedWords)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        if args.blackListedWords != "":
            instance.read_black_listed_words_file()           

        Timer.get_instance().step( "Read lncRNA-disease file..")
        instance.read_lncrna_disease_file( )
        
        Timer.get_instance().step( "Read protein-disease file..")
        instance.read_protein_disease_file( )
        
        Timer.get_instance().step( "Read lncRNA-protein enrichments file..")
        instance.read_lncrna_protein_file( )
        
        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

