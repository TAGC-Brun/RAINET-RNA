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
# Started 21-Oct-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to remove transcript redundancy by merging transcripts from same gene by choosing longest transcript only."
SCRIPT_NAME = "remove_transcript_redundancy.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read list of interacting RNAs (i.e. RNAs with interaction data) 
# 2) Read list of gene IDs, associated transcript interacting IDs and associated lengths 
# 2) For each gene, select a single transcript, the longest transcript, or the first to appear in case of draw
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Input file for this script can be created querying: SELECT transcriptID,geneID,transcriptLength FROM RNA
# 2) Decision to use flat files for input instead of a RAINET DB is that in this way we are sure to be using the right input files
#===============================================================================

# #
def run( list_interacting_rnas, info_file, outputGenes):

    #===============================================================================
    # Read list of interacting RNAs
    #===============================================================================

    listInteractingRnas = set()
    with open( list_interacting_rnas) as inFile:
        for line in inFile:
            line = line.strip()
            if "ensembID" in line:
                continue

            listInteractingRnas.add( line)

    print "Read %s interacting RNAs." % len( listInteractingRnas)

    #===============================================================================
    # Read geneID, txID, lengths
    #===============================================================================

    # Example:
    # "ENST00000002165","ENSG00000001036","2356"
    # "ENST00000002501","ENSG00000003249","2079"

    #     SQLManager.get_instance().set_DBpath( DBPATH)
    #     sql_session = SQLManager.get_instance().get_session()

    genes = {} # key -> geneID, value -> list of transcriptID + length

    nlines = 0
    with open( info_file) as inFile:
        for line in inFile:
            nlines += 1
            line = line.strip()
            spl = line.split(",")
            geneID = spl[1].replace('"','')
            transcriptID = spl[0].replace('"','')
            length = spl[2].replace('"','')

            # keep only genes/transcripts that have interaction information
            if transcriptID not in listInteractingRnas:
                continue
            
            if geneID not in genes:
                genes[ geneID] = []

            genes[ geneID].append( transcriptID + "|" + length)


    print "Number lines read: %s" % ( str(nlines))
    print "Number of genes with interacting transcripts: %s" % len( genes)

    if outputGenes:
        outFile = open("list_non_redundant_genes.txt" , "w")        
    else:
        outFile = open("list_non_redundant_transcripts.txt" , "w")

    count = 0    
    for geneID in genes:
        
        if outputGenes:
            outFile.write( geneID + "\n")
            count += 1
        else:        
            maxLength = 0
            maxTranscript = ""
            
            for tx in genes[ geneID]:
                spl = tx.split("|")
                length = int( spl[1])
                if length > maxLength:
                    maxLength = length
                    maxTranscript = spl[0]
    
            outFile.write( maxTranscript + "\n")
            count += 1
    
            if maxLength > 1200:
                print "Warning: length above 1200 nt"
                print maxLength, maxTranscript, geneID

    assert( count == len( genes))

    outFile.close()


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
        parser.add_argument('listOfInteractingRNAs', metavar='listOfInteractingRNAs', type=str,
                             help='List of interacting RNAs (RNAs with interaction data) from a catRAPID omics/all vs all file. Only these transcripts will be considered.')
        parser.add_argument('infoFile', metavar='infoFile', type=str,
                             help='Dump of RAINET DB RNA table. Tx ID, Gene ID, Tx length. e.g.: "ENST00000000233","ENSG00000004059","1103"')
        parser.add_argument('--outputGenes', metavar='outputGenes', type=int, default = 0,
                             help='Instead of outputting the longest transcript for a gene, output the gene ids. Default = 0 (OFF)')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        run( args.listOfInteractingRNAs, args.infoFile, args.outputGenes)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

