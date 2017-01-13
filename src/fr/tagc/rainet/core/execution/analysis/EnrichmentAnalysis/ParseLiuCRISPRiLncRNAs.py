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
# Started 13-January-2017 
# Diogo Ribeiro
DESC_COMMENT = "Script to read data from Liu2016 paper."
SCRIPT_NAME = "ParseLiuCRISPRIiLncRNAs.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read Table S1, containing correspondence between Liu2016 IDs and Ensembl transcript IDs
# 2) Read Table S10, containing information including effect of lncRNA on cell growth
# 3) Return list of Ensembl IDs that affect cell growth
#===============================================================================

#===============================================================================
# Processing notes:
# 1) processing both transcript IDs and gene IDs for sake of completion
#===============================================================================

class ParseLiuCRISPRIiLncRNAs(object):
    
    OUTPUT_FILE_TRANSCRIPTS = "/Liu2016_growth_affecting_lncRNAs_transcriptID.txt"
    OUTPUT_FILE_GENES = "/Liu2016_growth_affecting_lncRNAs_geneID.txt"
    
    def __init__(self, rainetDB, tableS1, tableS10, outputFolder):

        self.rainetDB = rainetDB
        self.tableS1 = tableS1
        self.tableS10 = tableS10
        self.outputFolder = outputFolder
        
        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)


    # #
    # Get correspondence between different RNA IDs using RainetDB
    def read_rainet_db(self):

        rnaCrossReference = {} # key -> cross reference ID, val -> ensembl Gene ID
 
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
            rnaCrossReference[ geneName].add( str( geneID) )
            # add external transcript name
            rnaCrossReference[ txName].add( str( geneID) )            
            # add gene-to-transcript correspondence
            rnaCrossReference[ geneID].add( str( geneID) )
            # add transcript entry
            rnaCrossReference[ ensemblID].add( str( geneID) )
  
        print "read_rainet_db: Read RNA table."
        print "read %s entries" % len( query)
 
        self.rnaCrossReference = rnaCrossReference


    # #
    # Read Table S1 and map Liu2016 IDs to Ensembl transcript IDs
    def read_tableS1(self):

        #===============================================================================
        # Read table S1
        #===============================================================================
        # Example format
        # Gene ID TSS ID  Transcript ID   gene name       chromosome      strand  TSS source      primary TSS 5prime      primary TSS 3prime      secondary TSS 5prime    secondary TSS 3prime
        # LH00001 TSS5    CUFF.182.5,CUFF.182.7   RP4-669L17.10   chr1    +       Annotation      319165  319165  319165  319165
        # LH00002 TSS11   ENST00000419160 RP4-669L17.10   chr1    +       Annotation      322672  322672  322672  322672

        # Note: they have data at the transcript level, and their experiments were performed at the transcript level.
        # they have different types of transcript IDs since they use different lncRNA annotation datasets in their study. These do not map to Ensembl transcript IDs.
        # they specifically used Ensembl v75 as one of their datasets, an old version containing only 15.574 transcripts from Ensembl.
        # therefore, we can have two approches:
        #     - Approach 1: map only transcript IDs from Ensembl (ENST*). Ignore all other transcriptIDs.
        #     - Approach 2: map all possible Gene names to ensembl Gene IDs (therefore making analysis at the gene level). This will include much more data.

        
        transcriptDict = {} # key -> Liu ID, value -> set of ensembl Transcript IDs
        ensemblSet = set()
       
        geneDict = {} # key -> Liu ID, value -> set of ensembl Gene IDs
        genesNotFound = set()
        genesFound = set()
        ensemblGeneSet = set()
        nLines = 0
 
        with open( self.tableS1, "r") as inFile:
            header = inFile.readline()
            for line in inFile:
                spl = line.strip().split( "\t")

                nLines += 1
                
                liuID = spl[0]
                # an entry can have several transcript IDs
                transcriptsIDs = spl[2].split( ",")
                # an entry can have several gene IDs
                geneIDs = spl[3].split( ",")

                ## Approach 1: direct transcript ENST correspondence                
                for tx in transcriptsIDs:
                    if tx.startswith("ENST"):
                        if liuID not in transcriptDict:
                            transcriptDict[ liuID] = set()
                        else:
                            # does not happen
                            pass

                        transcriptDict[ liuID].add( tx)
                        ensemblSet.add( tx)

                ## Approach 2: mapping gene name to ENSG
                for gene in geneIDs:
                    if gene in self.rnaCrossReference:
                        if liuID not in geneDict:
                            geneDict[ liuID] = set()
                        else:
                            # normal since there is repetition on the input file which is transcript-level
                            pass

                        # mapping can be to several genes, e.g. alternative chromosome assemblies etc
                        mappedGenes = self.rnaCrossReference[ gene]
                        for mapped in mappedGenes:
                            geneDict[ liuID].add( mapped)
                            ensemblGeneSet.add( mapped)

                        genesFound.add( gene)
                            
                    else:
                        genesNotFound.add( gene)
                        
        print "read_tableS1: Read %s lines." % nLines

        print "read_tableS1: %s Liu transcript IDs with Ensembl transcript ID correspondence." % len( transcriptDict)
        print "read_tableS1: %s Ensembl transcript IDs found." % len( ensemblSet)
        
        #grep -o ENST aah7111-TableS1.txt | wc -l 3623

        print "read_tableS1: %s Liu genes not found." % len( genesNotFound)
        print "read_tableS1: %s Liu genes found." % len( genesFound)
        print "read_tableS1: %s Liu gene IDs mapped to Ensembl gene IDs." % len( geneDict)
        print "read_tableS1: %s Ensembl gene IDs found." % len( ensemblGeneSet)

        self.transcriptDict = transcriptDict
        self.geneDict = geneDict


    # #
    # Read Table S10 and get list of lncRNA transcripts and genes affecting growth.
    # Writes output files
    def read_tableS10(self):

        #===============================================================================
        # Read table S10
        #===============================================================================
        # Example format
        # cell    gene    hit or non-hit  log2(FPKM)      TSS-pc distance Locus-locus distance    Transcript Length       Number of exons Is intergenic   Is antisense    Near Vista Enhancer     Near FANTOM enhancer
        #     Near Hnisz enhancer     Near Hnisz Super Enhancer       Near Cancer Associated SNP      Within Pol2 Loop        Within CTCF loop        Has mouse ortholog      Locus is amplified      Locus is heterozygous deleted   Locus is homozygous deleted     Locus is normal
        # K562    LH00001 0       2.9611936042    48473   43833   3019.5  2.5     1       0       0       0       0       0       0       0       0       0       0       0       0       1
        # K562    LH00002 0       -0.1795195586   44966   42685   609     2       1       0       0       0       0       0       0       0       0       0       0       0       0       1

        # Note: the column 3 (1-based, "hit or non-hit") corresponds to the association of this transcript to cell growth.
        # Note: there is data for several cell-types. We simply want to get the largest possible list of transcripts/genes, 
        # therefore being associated to growth in a single cell-type is enough.

        #===============================================================================
        # Output files
        #===============================================================================
        # simple lists of transcripts / genes
        outFile1 = open( self.outputFolder + ParseLiuCRISPRIiLncRNAs.OUTPUT_FILE_TRANSCRIPTS, "w")
        outFile2 = open( self.outputFolder + ParseLiuCRISPRIiLncRNAs.OUTPUT_FILE_GENES, "w")      

        nLines = 0

        with open( self.tableS10, "r") as inFile:
            header = inFile.readline()

            for line in inFile:
                spl = line.strip().split( "\t")
                
                nLines += 1
                
                liuGeneID = spl[1]
                hit = int( spl[2])

                # if transcript affecting cell growth   
                if hit:
                    
                    # check if there are ensembl transcripts mapped
                    if liuGeneID in self.transcriptDict:
                        transcripts = self.transcriptDict[ liuGeneID]
                        for tx in transcripts:
                            outFile1.write( tx + "\n")

                    # check if there are ensembl genes mapped                        
                    if liuGeneID in self.geneDict:
                        genes = self.geneDict[ liuGeneID]
                        for ge in genes:
                            outFile2.write( ge + "\n")

        print "read_tableS10: read %s lines." % nLines

        outFile1.close()
        outFile2.close()


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
                             help='Rainet database to be used for ID mapping.')
        parser.add_argument('tableS1', metavar='tableS1', type=str,
                             help='Liu2016 table S1, containing correspondence between their IDs and Ensembl IDs.')
        parser.add_argument('tableS10', metavar='tableS10', type=str,
                             help='Liu2016 table S10, containing information including effect of lncRNA on cell growth.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        parseLncrnadisease = ParseLiuCRISPRIiLncRNAs( args.rainetDB, args.tableS1, args.tableS10, args.outputFolder)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        Timer.get_instance().step( "Read Rainet DB..")
        parseLncrnadisease.read_rainet_db( )

        Timer.get_instance().step( "Read table S1..")
        parseLncrnadisease.read_tableS1( )

        Timer.get_instance().step( "Read table S10..")            
        parseLncrnadisease.read_tableS10( )


        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

