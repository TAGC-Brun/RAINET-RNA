
import sys
import os
import argparse

# import numpy as np
# import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

from Enrichment2Protein import Enrichment2Protein
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.ProteinRNAInteractionCatRAPID import ProteinRNAInteractionCatRAPID

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

    OUTPUT_FILE_WORD_MATCH = "/lncRNA_protein_disease_descriptions_word_match.txt"

    def __init__(self, rainetDB, lncRNADiseaseFile, proteinDiseaseFile, enrichmentData, outputFolder, minWordSize, blackListedWords, complexDatasets):

        self.rainetDB = rainetDB
        self.lncRNADiseaseFile = lncRNADiseaseFile
        self.proteinDiseaseFile = proteinDiseaseFile
        self.enrichmentData = enrichmentData
        self.outputFolder = outputFolder
        self.minWordSize = minWordSize
        self.blackListedWords = blackListedWords
        self.complexDatasets = complexDatasets

        # stores words to be black listed when performing word matching
        self.blackListedWordsSet = set()

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()


    # #
    # Get geneID-transcriptID correspondence and interactions for wanted transcripts
    def read_rainet_db(self):

        #===============================================================================
        # Query the RNA table to get transcriptID-geneID correspondence
        #===============================================================================

        geneTranscriptDict = {} # key -> transcriptID, value -> geneID

        for transcript in self.transcriptDiseaseDict:
            
            item = '"' + transcript + '"'

#             queryText = "query( RNA.geneID ).filter( RNA.transcriptID == %s).all()" % item
            queryText = "query( RNA.externalGeneName ).filter( RNA.transcriptID == %s).all()" % item
            
            res = eval('self.sql_session.' + queryText)

            geneTranscriptDict[ transcript] = str( res[0][0])
            

        #===============================================================================
        # Query the ProteinRNAInteractionCatRAPID to get interactions for wanted transcripts
        #===============================================================================

        transcriptInteractionDict = {} # key -> transcriptID, value -> set of interacting proteins

        for transcript in self.transcriptDiseaseDict:
            
            item = '"' + transcript + '"'

            queryText = "query( ProteinRNAInteractionCatRAPID.proteinID ).filter( ProteinRNAInteractionCatRAPID.transcriptID == %s).all()" % item
            
            res = eval('self.sql_session.' + queryText)
            
            if transcript not in transcriptInteractionDict:
                transcriptInteractionDict[ transcript] = set()
            
            for protein in res:
                transcriptInteractionDict[ transcript].add( str( protein[0]) )

        self.transcriptInteractionDict = transcriptInteractionDict          
        self.geneTranscriptDict = geneTranscriptDict
     
            
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
    # Map enrichment results to proteins of the complex, using another class
    def enrichment_to_protein(self):
        
        enrichment2Protein = Enrichment2Protein( self.rainetDB, self.enrichmentData, self.outputFolder, self.complexDatasets)
        
        pairsToWrite, enrichmentsDict = enrichment2Protein.run()

        #enrichmentsDict is # Key -> enrichment tag, val -> enrichment text

        self.pairsToWrite = pairsToWrite
        self.enrichmentsDict = enrichmentsDict
        
        return pairsToWrite, enrichmentsDict


    # #
    # Loop transcript-protein enrichment correspondence and write output file with disease associations
    def process_lncrna_protein_map(self):

        print "read_lncrna_protein_file: avoiding word match of words below size %s." % self.minWordSize
        
        matchInfo = [] # all information about the matched disease and its components

        # line counters after writing header
        nLines = 1
        nLinesMatch = 1        
        for enrich in self.enrichmentsDict:
            line = self.enrichmentsDict[ enrich]
            spl = line.strip().split( "\t")
                
            # e.g. line: ENST00000320202 44a     4       5       25.0%   0       9.6e-04 1.7e-02 1       BioplexCluster  Q9BVW5,Q9BRX5,Q9UNS1,Q9BRT9,Q14691,Q9Y248
    
            transcriptID = spl[0]
            
            if transcriptID not in self.transcriptDiseaseDict:
                # if transcript has no associated disease, skip transcript
                continue
            
            complexID = spl[1]
            datasetID = spl[9]
            enrichmentTag = complexID + "|" + datasetID
            # this entry contains several proteins in the complex, we want to see if at least one of them as disease overlap
            proteins = spl[10].split(",")
            
            geneID = self.geneTranscriptDict[ transcriptID]
                                    
                                
            # Match transcript and protein disease for this enrichment
            for transcriptDisease in self.transcriptDiseaseDict[ transcriptID]:

                ## Word matching
                # put all text in lowercase for word matching
                txDiseaseLower = transcriptDisease.lower()
                # Split per word
                wordSoup = txDiseaseLower.split( " ")
                # Filter out small words
                wordSoup = [word for word in wordSoup if len( word) > self.minWordSize]
                
                for proteinID in proteins:
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
                                    
                                    # get complexInteractions overlapping with complex
                                    proteinsInteracting = self.transcriptInteractionDict[ transcriptID].intersection( set( proteins))                                    
                                    assert len( proteinsInteracting) > 0
                                    proteinsInteractingText = ",".join( proteinsInteracting)
                                    
                                    matchInfo.append( "%s\t%s\t%s\t'%s'\t%s\t%s\t%s\t%s\t%s\n" % ( transcriptID, geneID, proteinID, matchWord, transcriptDisease, proteinDisease, enrichmentTag, spl[10], proteinsInteractingText))
                                    nLinesMatch += 1
                                                        
                            nLines += 1


        #===============================================================================
        # Output file, word matching
        #===============================================================================

        # File with matched disease terms only
        outFile = open( self.outputFolder + CommonLncRNAProteinDisease.OUTPUT_FILE_WORD_MATCH, "w")        
        outFile.write( "transcriptID\tGeneID\tproteinID\twordsMatched\ttranscriptDisease\tproteinDisease\tcomplexID\tcomplexProteins\tcomplexInteractions\n")
        
        for match in matchInfo:
            outFile.write( match )

        outFile.close()
        
        print "read_lncrna_protein_file: wrote %s association lines." % nLines
        print "read_lncrna_protein_file: wrote %s association lines with a word match." % nLinesMatch
        
        return matchInfo
        

    # #
    # Run all functions is order
    def run(self):
        
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        if self.blackListedWords != "":
            self.read_black_listed_words_file()           

        Timer.get_instance().step( "Read lncRNA-disease file..")
        self.read_lncrna_disease_file( )

        Timer.get_instance().step( "Read RainetDB table..")
        self.read_rainet_db( )
        
        Timer.get_instance().step( "Read protein-disease file..")
        self.read_protein_disease_file( )

        Timer.get_instance().step( "Read enrichment file..")        
        self.enrichment_to_protein()

        Timer.get_instance().step( "Process lncrna protein disease and enrichment map..")
        self.process_lncrna_protein_map()


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
        parser.add_argument('rainetDB', metavar='rainetDB', type=str,
                             help='Rainet database to be used.')
        parser.add_argument('lncRNADiseaseFile', metavar='lncRNADiseaseFile', type=str,
                             help='TSV file with lncRNA-disease associations. One line per pair. E.g. ENST00000428939 colorectal cancer')
        parser.add_argument('proteinDiseaseFile', metavar='proteinDiseaseFile', type=str,
                             help='TSV file with protein-disease associations (e.g. from OMIM). One line per pair. E.g. O00623  266510  PEROXISOME BIOGENESIS DISORDER 3B; PBD3B |  | ')
        parser.add_argument('enrichmentData', metavar='enrichmentData', type=str,
                             help='Data from enrichment analysis. This file should be pre-filtered for the wanted enrichments. Last column needs to be dataset of origin (e.g. WanCluster, BioplexCluster, etc), including change in header with "dataset" column.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
        parser.add_argument('--minWordSize', metavar='minWordSize', type=int, default = 2,
                             help='For matching transcript disease to protein disease, try to match any words above size X.')
        parser.add_argument('--blackListedWords', metavar='blackListedWords', type=str, default = "",
                             help='Optional file with list of words (one per line) that should not be matched by themselves, unless in conjunction with other words.')
        parser.add_argument('--complexDatasets', metavar='complexDatasets', type=str, default = "WanCluster,BioplexCluster,CustomCluster,NetworkModule",
                             help='List of datasets to withdrawn data from RAINET DB. Comma-separated.')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        instance = CommonLncRNAProteinDisease( args.rainetDB, args.lncRNADiseaseFile, args.proteinDiseaseFile, args.enrichmentData, args.outputFolder, args.minWordSize, args.blackListedWords, args.complexDatasets)
            
        instance.run()
            
        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

