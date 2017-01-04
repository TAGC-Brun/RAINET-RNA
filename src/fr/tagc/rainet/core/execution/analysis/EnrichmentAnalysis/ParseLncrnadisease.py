import sys
import os
import argparse

# import numpy as np
# import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

# from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
# from fr.tagc.rainet.core.util.data.DataManager import DataManager

from fr.tagc.rainet.core.data.RNACrossReference import RNACrossReference
from fr.tagc.rainet.core.data.RNA import RNA

#===============================================================================
# Started 03-January-2017 
# Diogo Ribeiro
DESC_COMMENT = "Script to read and ID map data from lncrnadisease database."
SCRIPT_NAME = "ParseLncrnadisease.py"
# Based on ParseLnc2cancer.py
#===============================================================================

#===============================================================================
# General plan:
# 1) Read Rainet DB containing tables containing RNA ID mapping
# 2) Read lncrnadisease file
# 3) Return lncrnadisease data with extra column harboring Ensembl IDs
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Use of a greedy approach that tries to map any lncRNA name / symbol / synonym to any ENST*
# 2) Output a line per transcript match, even if what is matched is the gene and not the transcript
#===============================================================================

class ParseLncrnadisease(object):
    
    OUTPUT_FILE = "/lncrnadisease_association_ensemblID.txt"
    
    def __init__(self, rainetDB, lncrnadiseaseData, outputFolder):

        self.rainetDB = rainetDB
        self.lncrnadiseaseData = lncrnadiseaseData
        self.outputFolder = outputFolder
        
        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)

    # #
    # Get correspondence of complex to protein IDs from RainetDB
    def read_rainet_db(self):

        #===============================================================================
        # Query CrossReference table
        #===============================================================================
        
        # query table containing all annotation mappings of rna and crossreference
        query = self.sql_session.query( RNACrossReference.transcriptID, RNACrossReference.crossReferenceID ).all()

        rnaCrossReference = {} # key -> cross reference ID, val -> ensembl transcript ID

        # a crossReference ID can map to several ensembl transcript IDs
        
        setOfTranscriptIDs = set()
        for ensemblID, crossRef in query:
            if crossRef not in rnaCrossReference:
                rnaCrossReference[ crossRef] = []
                
            rnaCrossReference[ crossRef].append( str( ensemblID) )

            setOfTranscriptIDs.add( ensemblID)

        print "read_rainet_db: Read RNACrossReference table."
        print "read %s entries" % len( query)
        print "read %s cross references for %s transcripts" % (len( rnaCrossReference), len(setOfTranscriptIDs) )

        #===============================================================================
        # Query RNA table, which contains gene ID, external gene and transcript names 
        #===============================================================================
        query = self.sql_session.query( RNA.transcriptID, RNA.geneID, RNA.externalGeneName, RNA.externalTranscriptName ).all()
        # Note: an gene name points to several ensembl IDs        
        
        for ensemblID, geneID, geneName, txName in query:
            if geneName not in rnaCrossReference:
                rnaCrossReference[ geneName] = set()
            if txName not in rnaCrossReference:
                rnaCrossReference[ txName] = set()
            if geneID not in rnaCrossReference:
                rnaCrossReference[ geneID] = set()
            if ensemblID not in rnaCrossReference:
                rnaCrossReference[ ensemblID] = set()

            # add external gene name                
            rnaCrossReference[ geneName].add( str( ensemblID) )
            # add external transcript name
            rnaCrossReference[ txName].add( str( ensemblID) )            
            # add gene-to-transcript correspondence
            rnaCrossReference[ geneID].add( str( ensemblID) )
            # add transcript entry
            rnaCrossReference[ ensemblID].add( str( ensemblID) )

            setOfTranscriptIDs.add( ensemblID)

        print "read_rainet_db: Read RNA table."
        print "read %s entries" % len( query)
        print "read %s cross references for %s transcripts" % (len( rnaCrossReference), len(setOfTranscriptIDs) )

        self.rnaCrossReference = rnaCrossReference


    # #
    # Read lncrnadisease data and map IDs to Ensembl transcript IDs
    def read_lncrnadisease(self):
        
        
        #===============================================================================
        # Reac lncrnadisease file
        #===============================================================================
        # Example format
        # Counter LncRNA name     Disease name    Dysfunction type        Functional description  Chr     Start   End     Strand  Species Synonyms        Other ID        PubMed ID
        # 2       1B FGF-antisense transcripts    endometriosis   expression      Mihalich et al, reported patients with endometriosis show low expression of 1B FGF-antisense transcripts, which correlates with endometrial cell proliferation  N/A     N/A     N/A     N/A     Human   1B FGF-antisense transcripts    N/A     23781896

        # Note: the original file contains no header, but input is a modified file with header

        # Approach: They provide lncRNA name and synonyms. 
        # First we try to map the lncRNA name and if this is not found, look for the synonyms and other IDs
        # Note: both lncRNA name and synonym columns are always filled. The synonym column often has the same as lncRNAName
        # Note: lncrnadisease provide one line per lncRNA-disease association, therefore, several lines may contain the same lncRNA associated to different diseases
        # Note: lncrnadisease IDs most often refer to Genes, here we map to transcripts, which may be all the transcripts of each of those genes
        # Note: lncrnadisease has data on different species, we are only interested in human
              
        notFound = set()
        found = set()

        #TODO: here
        
        # counters for mapping provinience
        namCount = 0
        synCount = 0
        refCount = 0        

        #===============================================================================
        # Output file
        #===============================================================================
        # Same as input file, but with an extra column (the first column) with the mapped Ensembl transcript ID
        outFile = open( self.outputFolder + ParseLncrnadisease.OUTPUT_FILE, "w")
    

        nLines = 0       
        with open( self.lncrnadiseaseData, "r") as inFile:

            # write header of output file            
            header = inFile.readline()
            newHeader = header.split("\t")
            # note that the first column of original header is empty. that entry is replaced with the new column
            newList = ["transcriptID"]
            newList.extend( newHeader[1:])
            newHeader = "\t".join( newList)
            outFile.write( newHeader)
                        
            for line in inFile:
                spl = line.strip().split( "\t")
                
                nLines += 1
                
                # first column of original file is just an incrementing number
                lncRNAName = spl[1]
                synonyms = spl[10].split( ";")                                                                                                                                                                                                                                   
                otherID = spl[4]                                                   
                species = spl[9]

                # filter out entries that are not in human                
                if species != "Human":
                    continue

                mappedTranscripts = set()                                                                   

                # Correct cases such as ENSG00000135253.9 to be ENSG00000135253
                if lncRNAName.startswith("ENSG") and "." in lncRNAName:
                    lncRNAName = lncRNAName.split(".")[0]
                    
                # 1) Try to map lncRNA name (i.e. mostly external gene name)
                if lncRNAName in self.rnaCrossReference:
                    for tx in self.rnaCrossReference[ lncRNAName]:
                        mappedTranscripts.add( tx)
                    found.add( lncRNAName)
                    namCount += 1

                else:
                    # 2) Try to map synonym lncRNA names / IDs

                    foundSyn = 0 # whether any synonym is found
                
                    for syn in synonyms:
                        syn = syn.strip() # they put a space after the ";"

                        if syn in self.rnaCrossReference:

                            for tx in self.rnaCrossReference[ syn]:
                                mappedTranscripts.add( tx)
                            
                            foundSyn = 1
                            found.add( lncRNAName)
                            synCount += 1
                            break

                    if foundSyn == 0:
                        # 3) Try to map Other IDs

                        if otherID in self.rnaCrossReference:                                    
                            for tx in self.rnaCrossReference[ otherID]:
                                mappedTranscripts.add( tx)
                            found.add( lncRNAName)
                            refCount += 1
                        else:
                            # if no mapping is found
                            notFound.add( lncRNAName)

                # write to output file
                if len( mappedTranscripts) > 0:
                    for transcript in mappedTranscripts:
                        newLine = [transcript]
                        newLine.extend( spl)
                        newText = "\t".join( newLine)
                        outFile.write( newText + "\n")

        outFile.close()

        print namCount, synCount, refCount

        print "read_lncrnadisease: read %s lines: " % nLines

        print "read_lncrnadisease: IDs found: %s" % len( found)

        print "read_lncrnadisease: IDs not found: %s" % len( notFound)
        print "List of transcripts not found: "
        print notFound
                        


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
        parser.add_argument('lncrnadiseaseData', metavar='lncrnadiseaseData', type=str,
                             help='Data from lncrnadisease database. This file must have a header.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        parseLncrnadisease = ParseLncrnadisease( args.rainetDB, args.lncrnadiseaseData, args.outputFolder)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        Timer.get_instance().step( "Read RainetDB table..")
        parseLncrnadisease.read_rainet_db( )

        Timer.get_instance().step( "Read lncrnadisease file..")            
        parseLncrnadisease.read_lncrnadisease( )


        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

