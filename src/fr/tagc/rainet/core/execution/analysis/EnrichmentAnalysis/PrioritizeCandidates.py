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

from sqlalchemy.inspection import inspect    

from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinInteraction import ProteinInteraction

from fr.tagc.rainet.core.data.BioplexCluster import BioplexCluster
from fr.tagc.rainet.core.data.ProteinBioplexAnnotation import ProteinBioplexAnnotation
from fr.tagc.rainet.core.data.WanCluster import WanCluster
from fr.tagc.rainet.core.data.ProteinWanAnnotation import ProteinWanAnnotation
from fr.tagc.rainet.core.data.CorumCluster import CorumCluster
from fr.tagc.rainet.core.data.ProteinCorumAnnotation import ProteinCorumAnnotation
from fr.tagc.rainet.core.data.CustomCluster import CustomCluster
from fr.tagc.rainet.core.data.ProteinCustomAnnotation import ProteinCustomAnnotation
from fr.tagc.rainet.core.data.NetworkModule import NetworkModule
from fr.tagc.rainet.core.data.NetworkModuleAnnotation import NetworkModuleAnnotation
from fr.tagc.rainet.core.data.ProteinNetworkModule import ProteinNetworkModule

#===============================================================================
# Started 23-November-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to prioritize transcript-protein interactions to be validated experimentally. To be run after enrichment analysis and result parsing."
SCRIPT_NAME = "prioritize_candidates.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read Rainet DB containing table used for enrichment analysis
# 2) Read list of transcript-complex pairs we want to look for
# 3) Read transcript-protein (experimental) interaction file
# 4) Return enriched transcript-complexes with known experimental interactions
#===============================================================================

#===============================================================================
# Processing notes:
# 1)
# 2)
#===============================================================================

class PrioritizeCandidates(object):
    
    # correspondance between the base table and associated data
    ANNOTATION_TABLES_DICT = {"NetworkModule" : "ProteinNetworkModule",
                              "ReactomePathway" : "ProteinReactomeAnnotation",
                              "KEGGPathway" : "ProteinKEGGAnnotation",
                              "BioplexCluster" : "ProteinBioplexAnnotation",
                              "WanCluster" : "ProteinWanAnnotation",
                              "CorumCluster" : "ProteinCorumAnnotation",
                              "CustomCluster" : "ProteinCustomAnnotation"}
       
    OUTPUT_INTERACTION_FILE = "enrichment_with_interactions.txt"
       
    
    def __init__(self, rainetDB, complexTable, listTxComplexPairs, listTxProteinInteractions, topTxComplexPairs):

        self.rainetDB = rainetDB
        self.complexTable = complexTable
        self.listTxComplexPairs = listTxComplexPairs
        self.listTxProteinInteractions = listTxProteinInteractions
        self.topTxComplexPairs = topTxComplexPairs
        
        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()

#         # make output folder
#         if not os.path.exists( self.outputFolder):
#             os.mkdir( self.outputFolder)



    # Get correspondence of complex to protein IDs from RainetDB
    def read_rainet_db(self):
        
        tableNameBase = self.complexTable
        tableNameAnnotation = PrioritizeCandidates.ANNOTATION_TABLES_DICT[ tableNameBase]

        # Get table primary key names 
        primaryKeys = eval("inspect( " + tableNameAnnotation + ").primary_key") 
        # e.g. change from ProteinKEGGAnnotation.keggPathway_id to keggPathway_id
        primaryKeys = [ str(pk).replace(tableNameAnnotation, "")[ 1:] for pk in primaryKeys]
         
        # query table containing all annotation mappings of pathway-protein
        proteinAnnotations = eval("self.sql_session.query( " + tableNameAnnotation + " ).all()")
  
        #===================================================================   
        # Process annotations
        #=================================================================== 
         
        complexAnnotDict = {}  # key -> complex id, value -> list of proteins IDs
 
        for annot in proteinAnnotations:
            complexID = str(eval("annot." + primaryKeys[ 0]))
            protID = str(eval("annot." + primaryKeys[ 1]))

            if complexID not in complexAnnotDict:
                complexAnnotDict[ complexID] = set()

            complexAnnotDict[ complexID].add( protID)

        Logger.get_instance().info( "read_rainet_db: read %s annotations" % len( complexAnnotDict ) )

        self.complexAnnotDict = complexAnnotDict

    # #
    # Read file with list of transcript and complex pairs (enriched)
    def read_tx_complex_file( self):
        
        # Example format
        #
        # transcriptID    annotID transcript_enrichments  annot_enrichments
        # ENST00000518553 46      4       1
        # ENST00000523556 31      4       2
        # ENST00000518553 31      4       2
        # ENST00000523556 25      4       3

        transcriptComplexEnrichments = {} # key -> transcript ID, value -> complex ID # for top X tx-complex pairs

        count = 0
        with open( self.listTxComplexPairs, "r") as inFile:
            header = inFile.readline()
            for line in inFile:
                count += 1
                
                spl = line.strip().split("\t")
                
                txID = spl[ 0]
                complexID = spl[ 1]
                
                if txID not in transcriptComplexEnrichments:
                    transcriptComplexEnrichments[ txID] = set()
                transcriptComplexEnrichments[ txID].add( complexID)
    
                # skip tx-complex pairs that are over the defined top (if defined)
                if self.topTxComplexPairs != -1 and count >= self.topTxComplexPairs:
                    break

        Logger.get_instance().info( "read_tx_complex_file: read %s top lines" % count )
        Logger.get_instance().info( "read_tx_complex_file: read %s enriched transcripts" % len( transcriptComplexEnrichments ) )
    
        self.transcriptComplexEnrichments = transcriptComplexEnrichments
   
   
    # Read protein-RNA interactions file with experimental interactions
    def read_protein_rna_interactions_file(self):

        # Example format
        # ENST00000600002|Q9UKA9  -1.7    1       162
        # ENST00000602495|Q92804  -3.22   1       34000

        transcriptExperimentalInteractions = {} # key -> transcript ID, value -> set of interacting proteins

        count = 0
        with open( self.listTxProteinInteractions, "r") as inFile:
            for line in inFile:
                spl = line.strip().split("\t")
                
                transcriptID, proteinID = spl[0].split( "|")
                
                if transcriptID not in transcriptExperimentalInteractions:
                    transcriptExperimentalInteractions[ transcriptID] = set()
                transcriptExperimentalInteractions[ transcriptID].add( proteinID)

                count += 1
        

        Logger.get_instance().info( "read_protein_rna_interactions_file: read %s interactions" % count )
        
        Logger.get_instance().info( "read_protein_rna_interactions_file: %s transcripts with interactions" % len( transcriptExperimentalInteractions) )

        self.transcriptExperimentalInteractions = transcriptExperimentalInteractions


    # Function to run over transcript-complex enrichments and return the ones that have at least one experimental interaction
    def produce_output(self):

        # Output file:
        # transcriptID\tComplexID\tListProteinsWithExperimentalInteraction

        outFile = open( PrioritizeCandidates.OUTPUT_INTERACTION_FILE, "w")

        # loop over pairs of transcripts and their complex enrichments        
        for txID in self.transcriptComplexEnrichments:

            # is there any experimental interaction for this transcript
            if txID not in self.transcriptExperimentalInteractions:
                continue # not a priority candidate

            # experimentally interacting proteins
            interactingProteins = self.transcriptExperimentalInteractions[ txID]

            # for each enriched complex of this transcript
            for complexID in self.transcriptComplexEnrichments[ txID]:
                # list of proteins in current complex
                complexedProteins = self.complexAnnotDict[ complexID]   
                
                # overlap between proteins in complex and experimental interactions
                intersection = interactingProteins.intersection( complexedProteins)
                
                if len( intersection) > 0:
                    outFile.write( "%s\t%s\t%s\n" % ( txID, complexID, ",".join( intersection)) )  
    

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
        parser.add_argument('complexTable', metavar='complexTable', type=str,
                             help='Table within rainetDB to be used.')
        parser.add_argument('listTxComplexPairs', metavar='listTxComplexPairs', type=str,
                             help='List of transcript-complex pairs to be searched. E.g. enrichment_specificity_rank.tsv. File has header.')
        parser.add_argument('listTxProteinInteractions', metavar='listTxProteinInteractions', type=str,
                             help='List of (experimental) transcript-protein interactions, e.g. from NPInter, StarBase. E.g. ENST00000619449|Q96PU8.')        
        # optional args
        parser.add_argument('--topTxComplexPairs', metavar='topTxComplexPairs', type=int, default = -1,
                             help='Use only top X transcript-complex pairs while searching for candidates. (Default = -1, i.e. OFF).')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        prioritizeCandidates = PrioritizeCandidates( args.rainetDB, args.complexTable, args.listTxComplexPairs, args.listTxProteinInteractions, args.topTxComplexPairs)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        Timer.get_instance().step( "Read RainetDB table..")
        prioritizeCandidates.read_rainet_db( )

        Timer.get_instance().step( "Read transcript complex pairs file..")            
        prioritizeCandidates.read_tx_complex_file( )

        Timer.get_instance().step( "Read protein-RNA interactions file..")            
        prioritizeCandidates.read_protein_rna_interactions_file( )

        Timer.get_instance().step( "Produce output file..")            
        prioritizeCandidates.produce_output()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

