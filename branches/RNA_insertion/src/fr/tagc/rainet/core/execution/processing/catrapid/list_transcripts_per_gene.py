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
# Started 07-Mars-2017
# Diogo Ribeiro
DESC_COMMENT = "Script to produce list of transcripts, given a list of genes."
SCRIPT_NAME = "list_transcripts_per_gene.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read list of gene IDs  
# 2) Read gene-transcript correspondence
# 2) Output list of transcript IDs
#===============================================================================

# #
def run( gene_list, info_file):

    #===============================================================================
    # Read list of interacting RNAs
    #===============================================================================

    listGenes = set()
    with open( gene_list) as inFile:
        for line in inFile:
            line = line.strip()
            if "ensembID" in line:
                continue

            listGenes.add( line)

    print "Read %s genes." % len( listGenes)

    #===============================================================================
    # Read geneID, txID, lengths
    #===============================================================================

    # Example:
    # "ENST00000002165","ENSG00000001036","2356"
    # "ENST00000002501","ENSG00000003249","2079"

    #     SQLManager.get_instance().set_DBpath( DBPATH)
    #     sql_session = SQLManager.get_instance().get_session()

    outFile = open("list_transcripts_per_gene.txt","w")

    written = 0
    nlines = 0
    with open( info_file) as inFile:
        for line in inFile:
            nlines += 1
            line = line.strip()
            spl = line.split(",")
            geneID = spl[1].replace('"','')
            transcriptID = spl[0].replace('"','')

            if geneID not in listGenes:
                continue
            else:
                outFile.write( transcriptID + "\n")
                written+=1
            

    print "Number infoFile lines read: %s" % ( str(nlines))
    print "Wrote %s transcript IDs to %s" % ( str(written), outFile.name)


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
        parser.add_argument('geneList', metavar='geneList', type=str,
                             help='List of Ensembl gene IDs for which we want to know their transcripts.')
        parser.add_argument('infoFile', metavar='infoFile', type=str,
                             help='Dump of RAINET DB RNA table. Tx ID, Gene ID, Tx length. e.g.: "ENST00000000233","ENSG00000004059","1103"')
        
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        run( args.geneList, args.infoFile)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

