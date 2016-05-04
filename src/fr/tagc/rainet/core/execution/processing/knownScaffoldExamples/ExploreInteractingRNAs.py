import sys
import os
import argparse
import glob
import numpy as np
import random 
import pandas as pd
from scipy import stats

from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil

from sqlalchemy import or_,and_
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinInteraction import ProteinInteraction
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.RNACrossReference import RNACrossReference
from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.RNATissueExpression import RNATissueExpression
from fr.tagc.rainet.core.data.ProteinRNAInteractionCatRAPID import ProteinRNAInteractionCatRAPID
from fr.tagc.rainet.core.util.data.DataManager import DataManager

#===============================================================================
# Started 29-Apr-2016 
# Diogo Ribeiro
# Script to explore characteristics of a selected list of RNAs (Ensembl IDs)
#===============================================================================

#===============================================================================
# General plan:
#===============================================================================
#===============================================================================
# Processing notes:
#===============================================================================

class ExploreInteractingRNAs( object ):

    RNA_ALL_KW = "allRNAs"

    def __init__(self, rainetDB, transcriptList, interactingTranscriptList, outputFolder):

        self.rainetDB = rainetDB
        self.transcriptList = transcriptList
        self.interactingTranscriptList = interactingTranscriptList
        self.outputFolder = outputFolder
    
        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()


    # #
    # Read list of transcript ensembl IDs, one per line    
    def read_transcript_list(self, in_file):
    
        with open( in_file) as f:
            return [ line.strip() for line in f]
    
    # #
    #
    def query_database(self, transcript_list, interacting_transcript_list):

        # query all RNAs in database
        queryText = "query( RNA).all()"
        DataManager.get_instance().perform_query( ExploreInteractingRNAs.RNA_ALL_KW, queryText)
                
        # Convert format to dictionary with transcriptID as key
        DataManager.get_instance().query_to_object_dict( ExploreInteractingRNAs.RNA_ALL_KW, "transcriptID")
        allRNAs = DataManager.get_instance().get_data( ExploreInteractingRNAs.RNA_ALL_KW)

        # check list of RNAs
        notFound = set()
        transcriptList = set()
        for listedRNA in transcript_list:
            if listedRNA not in allRNAs:
                # Note: this is because catRAPID interactions file used different Ensembl version than the one in RAINET DB
                # Logger.get_instance().warning( "ExploreInteractingRNAs.query_database: listed RNA not found in database. %s" % listedRNA)
                # raise RainetException( "ExploreInteractingRNAs.query_database: listed RNA not found in database.", listedRNA)
                notFound.add( listedRNA)
            else:
                transcriptList.add( listedRNA)
        print "Listed RNAs not found in database:",len(notFound)
        print "Listed RNAs found in database:",len(transcriptList)
        
        notFoundInteracting = set()
        interactingTranscriptList = set()
        for listedRNA in interacting_transcript_list:
            if listedRNA not in allRNAs:
                # Note: this is because catRAPID interactions file used different Ensembl version than the one in RAINET DB
                # Logger.get_instance().warning( "ExploreInteractingRNAs.query_database: listed RNA not found in database. %s" % listedRNA)
                # raise RainetException( "ExploreInteractingRNAs.query_database: listed RNA not found in database.", listedRNA)
                notFoundInteracting.add( listedRNA)
            else:
                interactingTranscriptList.add( listedRNA)
        print "Listed interacting RNAs not found in database:",len(notFoundInteracting)
        print "Listed interacting RNAs found in database:",len(interactingTranscriptList)


        outFile = open(self.outputFolder+"/outFile.txt","w")
        outFile.write( "transcriptID\ttranscriptLength\ttranscriptBiotype\tinSet\n" )
        
        file2 = open(self.outputFolder+"/outFile2.txt","w")
        file3 = open(self.outputFolder+"/outFile3.txt","w")
        file3.write( "transcriptID\ttranscriptLength\ttranscriptBiotype\tinSet\n" )
        
        for rnaID in interactingTranscriptList:
            rna = allRNAs[ rnaID]

            lastCol = 0
            if rnaID in transcriptList:
                lastCol = 1
  
            outFile.write( "%s\t%s\t%s\t%s\n" % (rna.transcriptID, rna.transcriptLength, rna.transcriptBiotype, lastCol) )
            
            if rna.transcriptBiotype == "lincRNA":
                file2.write(rna.transcriptID+"\n")
                file3.write( "%s\t%s\t%s\t%s\n" % (rna.transcriptID, rna.transcriptLength, rna.transcriptBiotype, lastCol) )
             
        outFile.close()
 

if __name__ == "__main__":
    
    try:
        # Create Logger instance by using the first log action.
        Logger.get_instance().info( "ExploreInteractingRNAs : Starting..." )

        #===============================================================================
        # Get input arguments, initialise class
        #===============================================================================
        parser = argparse.ArgumentParser(description='# Script to explore properties of a list of RNAs.') 

        # positional args
        parser.add_argument('rainetDB', metavar='rainetDB', type=str, help='Path to RAINET database to be used.')
        parser.add_argument('transcriptList', metavar='transcriptList', type=str, help='List of Ensembl IDs, one per line.')
        parser.add_argument('interactingTranscriptList', metavar='interactingTranscriptList', type=str, help='List of background Ensembl IDs, one per line.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
        
        #gets the arguments
        args = parser.parse_args( ) 

        # Initialise class
        run = ExploreInteractingRNAs( args.rainetDB, args.transcriptList, args.interactingTranscriptList, args.outputFolder)

        #===============================================================================
        # Run analysis / processing
        #===============================================================================
         
        # Start chrono
        Timer.get_instance().start_chrono()
 
        transcriptList = run.read_transcript_list( args.transcriptList)
        interactingTranscriptList = run.read_transcript_list( args.interactingTranscriptList)

        run.query_database( transcriptList, interactingTranscriptList)
        


    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of ExploreInteractingRNAs. Aborting :\n" + rainet.to_string())

    # Stop the chrono      
    Timer.get_instance().stop_chrono( "ExploreInteractingRNAs : Finished" )


