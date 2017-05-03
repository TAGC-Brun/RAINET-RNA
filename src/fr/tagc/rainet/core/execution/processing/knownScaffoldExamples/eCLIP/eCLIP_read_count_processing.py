import sys
import os
import argparse
import glob

# import numpy as np
import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.RNACrossReference import RNACrossReference

# from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil

#===============================================================================
# Started 03-May-2017
# Diogo Ribeiro
DESC_COMMENT = "Script to process eCLIP mapped files containing read counts."
SCRIPT_NAME = "eCLIP_read_count_processing.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read RAINET DB for RefSeq-Ensembl mapping, as well as Ensembl Tx ID to Ensembl Gene ID correspondence.
# 2) Read eCLIP mapped files, create uniprotAC-EnsemblTxID interactions, containing read counts.
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Original files were mapped to RefSeq transcripts, an arbitrary transcript from a gene. Here we map them to Ensembl tx IDs 
# 2)
#===============================================================================


# 0-based columns
TX_ID_COLUMN = 0
READ_COUNT_COLUMN = 4
OUTPUT_INTERACTIONS = "eclip_interactions_readcount.out"

# #
# Get correspondence between RefSeq and Ensembl transcript IDs
# Get correspondence between transcript IDs and gene IDs using rainetDB
def read_rainet_db( sqlSession):


    #===============================================================================
    # Query RNA table, which contains gene ID
    #===============================================================================
    query = sqlSession.query( RNA.transcriptID, RNA.geneID ).all()
    # Note: an gene name points to several ensembl IDs        

    txGeneDict = {} # key -> ensembl transcript ID, val -> ensembl Gene ID
    geneTxDict = {} # key -> ensembl Gene ID, val -> list of ensembl transcript IDs
     
    # Correspondence should be many-to-1 (a transcript can only have an associated gene, a gene can have many transcripts) 

    for transcriptID, geneID in query:
        if transcriptID in txGeneDict:
            raise RainetException("read_rainet_db: duplicate transcript ID ")
        
        txGeneDict[ transcriptID] = str( geneID)

        if geneID not in geneTxDict:
            geneTxDict[ geneID] = set()
        geneTxDict[ geneID].add( transcriptID)

    Logger.get_instance().info( "read_rainet_db : Number Ensembl transcripts read %s" % ( len( txGeneDict)) )
    Logger.get_instance().info( "read_rainet_db : Number Ensembl genes read %s" % ( len( geneTxDict)) )
 
      
    #===============================================================================
    # Query CrossReference table, which contains RefSeq IDs
    #===============================================================================
    # For ncRNA, IDs belong to refseq_ncrna database
    # Not all refSeq ncRNAs have a Ensembl counterpart
    # A refSeq ID can have several ensembl transcript IDs
    query = sqlSession.query( RNACrossReference.transcriptID, RNACrossReference.crossReferenceID ).filter(RNACrossReference.sourceDB == "refseq_ncrna").all()

    refseqEnsemblDict = {} # key -> RefSeqID, value -> list of ensembl transcript ID
    setEnsemblIDs = set()

    for ensemblID, refseqID in query:
        
        if refseqID not in refseqEnsemblDict:
            refseqEnsemblDict[ refseqID] = set()
        refseqEnsemblDict[ refseqID].add( ensemblID)

        setEnsemblIDs.add( ensemblID)

    Logger.get_instance().info( "read_rainet_db : Number RefSeq transcripts read %s" % ( len( refseqEnsemblDict)) )
    Logger.get_instance().info( "read_rainet_db : Number Ensembl transcripts corresponding to RefSeq %s" % ( len( setEnsemblIDs)) )

    return txGeneDict, geneTxDict, refseqEnsemblDict


# #
# Read mapped eclip files and write interaction dataset using ensembl transcript IDs
def read_mapped_eclip_folder( eclipFolder, txGeneDict, geneTxDict, refseqEnsemblDict):


    # Processing notes:
    # In a eclip mapped folder, there are one file per RBP.
    # Name of file is uniprotAC
    # Each line in each file is mapping to a transcript, example format:
    # NM_001143984 14 26  20 2568 109.206
    # NM_001143985 4  0 2 1122 75.1204
    # NM_001003805 3  0 1.5  -1 11.5778
    # NM_001003803 3  0 1.5  -1 11.5778
    # NM_001003800 22  0 11 6433 5.70603

    # RefSeqTxID read_count_rep1 read_count_rep2 (space) mean_read_count length_of_transcript expression_level_of_transcript_in_cell_line 

    # Output file
    # TranscriptID\tProteinID\tclipReads
    outFile = open( OUTPUT_INTERACTIONS, "w")
    outFile.write( "transcript\tprotein\tclipReads\n")

    filesToProcess = glob.glob( eclipFolder + "/*")

    linesWritten = 0

    for file in filesToProcess:
        protID = file.split("/")[-1]
        
        inFile = open( file, "r")
        for line in inFile:
            line = line.strip()
            spl = line.split(" ")
            refSeqID = spl[ TX_ID_COLUMN]
            
            try:
                meanReads = float( spl[ READ_COUNT_COLUMN])
            except:
                print line
                print spl

            # Get Ensembl transcript ID from RefSeq ID
            setOfEnsemblTx = set()
            if refSeqID in refseqEnsemblDict:
                for ensemblID in refseqEnsemblDict[ refSeqID]:
                    setOfEnsemblTx.add( ensemblID)
# testing
#                     geneID = txGeneDict[ ensemblID]
#                     for geneTx in geneTxDict[ geneID]:
#                         setOfEnsemblTx.add( geneTx)

            # write interactions
            if len( setOfEnsemblTx) > 0:
                for tx in setOfEnsemblTx:
                    outFile.write( "%s\t%s\t%s\n" % (tx, protID, meanReads) )
                    linesWritten += 1
        
    outFile.close()

    Logger.get_instance().info( "read_mapped_eclip_folder : Wrote %s interactions." % ( linesWritten) )

        

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
        parser.add_argument('eclipFolder', metavar='eclipFolder', type=str,
                             help='Folder with mapped eclip files to be processed.')
        parser.add_argument('rainetDB', metavar='rainetDB', type=str, help='Path to RAINET database to be used.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(args.rainetDB)
        sqlSession = SQLManager.get_instance().get_session()

        Timer.get_instance().step( "Read RAINET DB..")            
        txGeneDict, geneTxDict, refseqEnsemblDict = read_rainet_db( sqlSession)

        Timer.get_instance().step( "Read mapped eCLIP folder..")            
        read_mapped_eclip_folder( args.eclipFolder, txGeneDict, geneTxDict, refseqEnsemblDict)



        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

