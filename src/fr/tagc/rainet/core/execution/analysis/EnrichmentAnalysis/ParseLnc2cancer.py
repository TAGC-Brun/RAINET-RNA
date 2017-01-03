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
DESC_COMMENT = "Script to read and ID map data from lnc2cancer database."
SCRIPT_NAME = "ParseLnc2cancer.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read Rainet DB containing tables containing RNA ID mapping
# 2) Read lnc2cancer file
# 3) Return lnc2cancer data with extra column harboring Ensembl IDs
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Use of a greedy approach that tries to map any lncRNA name / symbol / synonym to any ENST*
# 2) Output a line per transcript match, even if what is matched is the gene and not the transcript
#===============================================================================

class ParseLnc2cancer(object):
    
       
    
    def __init__(self, rainetDB, lnc2cancerData, outputFolder):

        self.rainetDB = rainetDB
        self.lnc2cancerData = lnc2cancerData
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
            # add transcript entry as well as some lnc2cancer IDs have ENST*
            rnaCrossReference[ ensemblID].add( str( ensemblID) )

            setOfTranscriptIDs.add( ensemblID)

        print "read_rainet_db: Read RNA table."
        print "read %s entries" % len( query)
        print "read %s cross references for %s transcripts" % (len( rnaCrossReference), len(setOfTranscriptIDs) )

        self.rnaCrossReference = rnaCrossReference


    # #
    # Read lnc2cancer data and map IDs to Ensembl transcript IDs
    def read_lnc2cancer(self):
        
        
        #===============================================================================
        # Reac lnc2cancer file
        #===============================================================================
        # Example format
        # [internal ID] LncRNA name     Synonyms        Ensembl ID      Refseq ID       Position        Cancer name     ICD-0-3-T       ICD-0-3-M       Methods Sample (tissue/ cell line)      Expression pattern      Functional description  PubMed ID       Year    Title
        # 1       91H     AI747191; 91H   N/A     N/A     N/A     colorectal cancer       C19.9           qPCR, RNAi etc. CRC tissue, cell lines ( HCT8, HT29, HCT116, SW620 etc.)        up-regulated    91H was significantly overexpressed in cancerous tissue and CRC cell lines compared with adjacent normal tissue and a normal human intestinal epithelial cell line. Moreover, 91H overexpression was closely associated with distant metastasis and poor prognosis in patients with CRC, except for CNV of 91H. Multivariate analysis indicated that 91H expression was an independent prognostic indicator, as well as distant metastasis. 91H played an important role in the molecular etiology of CRC and might be regarded as a novel prognosis indicator in patients with CRC.    25058480        2014    Up-regulation of 91H promotes tumor metastasis and predicts poor prognosis for patients with colorectal cancer.        
        # 10    ABHD11-AS1    ABHD11-AS1; WBSCR26; LINC00035; NCRNA00035    ENSG00000225969    NR_026690     chr7:73735038-73736054    gastric cancer    C16        qPCR etc.    gastric cancer tissue    up-regulated    Results show that compared with adjacent nontumor tissues the expression level of ABHD11-AS1 in gastric cancer tissues was significantly increased.     24984296    2014    Increased expression of long noncoding RNA ABHD11-AS1 in gastric cancer and its clinical significance.

        # Approach: They provide lncRNA name and synonyms. 
        # First we try to map the lncRNA name and if this is not found, look for the synonyms, then Ensembl ID, then RefSeq ID
        # Note: both columns are always filled. The synonym column has at least the same as the lncRNAName
        # Note: lnc2cancer provide one line per lncRNA-disease association, therefore, several lines may contain the same lncRNA associated to different diseases
        # Note: lnc2cancer IDs most often refer to Genes, here we map to transcripts, which may be all the transcripts of each of those genes
              
        notFound = set()
        found = set()
        
        # counters for mapping provinience
        namCount = 0
        synCount = 0
        ensCount = 0
        refCount = 0        

        #===============================================================================
        # Output file
        #===============================================================================
        # Same as input file, but with an extra column (the first column) with the mapped Ensembl transcript ID
        outFile = open( self.outputFolder + "/lncRNA_cancer_association_ensemblID.txt", "w")
    
       
        with open( self.lnc2cancerData, "r") as inFile:

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

                # first column of original file is just an incrementing number
                lncRNAName = spl[1]
                synonyms = spl[2].split( ";")
                ensemblID = spl[3]
                refseqID = spl[4]

                mappedTranscripts = set()

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
                        # 3) Try to map Ensembl IDs

                        if ensemblID in self.rnaCrossReference:                                    
                            for tx in self.rnaCrossReference[ ensemblID]:
                                mappedTranscripts.add( tx)
                            found.add( lncRNAName)
                            ensCount += 1
                        else:
                            # 4) Try to map RefSeq IDs
                            if refseqID in self.rnaCrossReference:                                    
                                for tx in self.rnaCrossReference[ refseqID]:
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

        # print namCount, synCount, ensCount, refCount

        print "read_lnc2cancer: lnc2cancer IDs found: %s" % len( found)

        print "read_lnc2cancer: lnc2cancer IDs not found: %s" % len( notFound)
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
        parser.add_argument('lnc2cancerData', metavar='lnc2cancerData', type=str,
                             help='Data from lnc2cancer database.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        prioritizeCandidates = ParseLnc2cancer( args.rainetDB, args.lnc2cancerData, args.outputFolder)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        Timer.get_instance().step( "Read RainetDB table..")
        prioritizeCandidates.read_rainet_db( )

        Timer.get_instance().step( "Read lnc2cancer file..")            
        prioritizeCandidates.read_lnc2cancer( )


        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

