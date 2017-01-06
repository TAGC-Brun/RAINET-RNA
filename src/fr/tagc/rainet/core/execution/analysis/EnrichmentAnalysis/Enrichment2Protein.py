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

from fr.tagc.rainet.core.data.DBParameter import DBParameter
from fr.tagc.rainet.core.data.GeneOntology import GeneOntology
from fr.tagc.rainet.core.data.GeneSymbol import GeneSymbol
from fr.tagc.rainet.core.data.KEGGPathway import KEGGPathway
from fr.tagc.rainet.core.data.NetworkModuleAnnotation import NetworkModuleAnnotation
from fr.tagc.rainet.core.data.NetworkModule import NetworkModule
from fr.tagc.rainet.core.data.OCGPartitionAnalysis import OCGPartitionAnalysis
from fr.tagc.rainet.core.data.PartitionAnalysis import PartitionAnalysis
from fr.tagc.rainet.core.data.PPINetworkInteraction import PPINetworkInteraction
from fr.tagc.rainet.core.data.PPINetwork import PPINetwork
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.ProteinDomain import ProteinDomain
from fr.tagc.rainet.core.data.ProteinGOAnnotation import ProteinGOAnnotation
from fr.tagc.rainet.core.data.ProteinInteraction import ProteinInteraction
from fr.tagc.rainet.core.data.ProteinIsoform import ProteinIsoform
from fr.tagc.rainet.core.data.ProteinKEGGAnnotation import ProteinKEGGAnnotation
from fr.tagc.rainet.core.data.ProteinNetworkModule import ProteinNetworkModule
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinReactomeAnnotation import ProteinReactomeAnnotation
from fr.tagc.rainet.core.data.ReactomePathway import ReactomePathway
from fr.tagc.rainet.core.data.SynonymGeneSymbol import SynonymGeneSymbol
from fr.tagc.rainet.core.data.Gene import Gene
from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.MRNA import MRNA
from fr.tagc.rainet.core.data.LncRNA import LncRNA
from fr.tagc.rainet.core.data.OtherRNA import OtherRNA
from fr.tagc.rainet.core.data.RNACrossReference import RNACrossReference
from fr.tagc.rainet.core.data.ProteinRNAInteractionCatRAPID import ProteinRNAInteractionCatRAPID
from fr.tagc.rainet.core.data.RNATissueExpression import RNATissueExpression
from fr.tagc.rainet.core.data.Tissue import Tissue
from fr.tagc.rainet.core.data.TableStatus import TableStatus
from fr.tagc.rainet.core.data.InteractingProtein import InteractingProtein
from fr.tagc.rainet.core.data.InteractingRNA import InteractingRNA
from fr.tagc.rainet.core.data.BioplexCluster import BioplexCluster
from fr.tagc.rainet.core.data.ProteinBioplexAnnotation import ProteinBioplexAnnotation
from fr.tagc.rainet.core.data.WanCluster import WanCluster
from fr.tagc.rainet.core.data.ProteinWanAnnotation import ProteinWanAnnotation
from fr.tagc.rainet.core.data.CorumCluster import CorumCluster
from fr.tagc.rainet.core.data.ProteinCorumAnnotation import ProteinCorumAnnotation
from fr.tagc.rainet.core.data.CustomCluster import CustomCluster
from fr.tagc.rainet.core.data.ProteinCustomAnnotation import ProteinCustomAnnotation


#===============================================================================
# Started 04-January-2017
# Diogo Ribeiro
DESC_COMMENT = "Script to map from enrichment results to having proteins of the complex."
SCRIPT_NAME = "Enrichment2Protein.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read (filtered) file with enrichments (with extra column with dataset of origin)
# 2) Read RAINET database complex-protein mapping
# 3) Return same file but with list of proteins associated to that complex (regardless of interaction)
#===============================================================================

#===============================================================================
# Processing notes:
# 1) The words 'complex' and 'cluster' are used interchangeably
#===============================================================================

class Enrichment2Protein(object):
    
    # Constants
    
    OUTPUT_FILE_1 = "/enrichment_results_with_protein_ids.txt"
    OUTPUT_FILE_2 = "/transcript_ids_associated_protein_ids.txt"

    # correspondance between the base table and associated data
    ANNOTATION_TABLES_DICT = {"NetworkModule" : "ProteinNetworkModule",
                              "ReactomePathway" : "ProteinReactomeAnnotation",
                              "KEGGPathway" : "ProteinKEGGAnnotation",
                              "BioplexCluster" : "ProteinBioplexAnnotation",
                              "WanCluster" : "ProteinWanAnnotation",
                              "CorumCluster" : "ProteinCorumAnnotation",
                              "CustomCluster" : "ProteinCustomAnnotation"}
    
    def __init__(self, rainetDB, enrichmentData, outputFolder, complexDatasets):

        self.rainetDB = rainetDB
        self.enrichmentData = enrichmentData
        self.outputFolder = outputFolder
        self.complexDatasets = complexDatasets.split(",")
        
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
        # Query the several complex-related tables
        #===============================================================================

        complexData = {} # Key -> cluster dataset, val -> dict. Key -> cluster ID, val -> list of protein IDs of that cluster

        # Loop over the several datasets
        for dataset in self.complexDatasets:

            # Define wanted tables         
            tableNameBase = dataset
            tableNameAnnotation = Enrichment2Protein.ANNOTATION_TABLES_DICT[ tableNameBase]
             
            # Get table primary key names 
            primaryKeys = eval("inspect( " + tableNameAnnotation + ").primary_key") 
            # e.g. change from ProteinKEGGAnnotation.keggPathway_id to keggPathway_id
            primaryKeys = [ str(pk).replace(tableNameAnnotation, "")[ 1:] for pk in primaryKeys]
             
            # query table containing all pathway/module descriptions
            #baseAnnotations = eval("self.sql_session.query( " + tableNameBase + " ).all()")
            # query table containing all annotation mappings of pathway-protein
            proteinAnnotations = eval("self.sql_session.query( " + tableNameAnnotation + " ).all()")

            if dataset not in complexData:
                complexData[ dataset] = {}
        
            for annot in proteinAnnotations:
                clusterID = str(eval("annot." + primaryKeys[ 0]))
                protID = str(eval("annot." + primaryKeys[ 1]))

                if clusterID not in complexData[ dataset]:
                    complexData[ dataset][ clusterID] = set()

                complexData[ dataset][ clusterID].add( protID)

            print "read_rainet_db: read %s which has %s entries." % ( dataset, len( complexData[ dataset]))
        
        self.complexData = complexData


    # #
    # Read enrichment table, produce output file
    def read_enrichment_data(self):
        
        #===============================================================================
        # Read enrichment file
        #===============================================================================
        # Example format
        # transcriptID    annotID number_observed_interactions    number_possible_interactions    percent_interacting_proteins    warning pval    corrected_pval  sign_corrected  dataset
        # ENST00000320202 44a     4       5       25.0%   0       9.6e-04 1.7e-02 1       BioplexCluster        
        # ENST00000610143 149     7       13      19.3%   0       9.0e-04 1.4e-02 1       CustomCluster


        #===============================================================================
        # Output files
        #===============================================================================
        # 1) as input file, but with extra column with list of protein IDs of that complex.
        # 2) Simple list of transcriptID-proteinID associations, one line per transcriptID-proteinID pair

        outFile1 = open( self.outputFolder + Enrichment2Protein.OUTPUT_FILE_1, "w")
        outFile2 = open( self.outputFolder + Enrichment2Protein.OUTPUT_FILE_2, "w")

        nLines = 0

        pairsToWrite = set()

        with open( self.enrichmentData, "r") as inFile:

            # get header and rewrite it with new column
            header = inFile.readline()
            newHeader = header.strip() + "\t" + "protein_list" + "\n"
            outFile1.write( newHeader)
            nLines += 1

            for line in inFile:
                spl = line.strip().split( "\t")

                nLines += 1
                
                transcriptID = spl[0]
                annotID = spl[1]
                dataset = spl[-1]
                clusterSize = int( spl[3])
                                
                # retrieve list of proteins on that dataset and complex
                proteinList = self.complexData[ dataset][ annotID]

                assert clusterSize <= len( proteinList), "number of protein list needs to be equal or higher than number of proteins with interactions"

                # write to novel enrichment file
                proteinText = ",".join( proteinList)                
                spl.extend( [proteinText])                
                newLine = "\t".join( spl)
                
                outFile1.write( newLine + "\n")

                # write list of transcriptID and proteinIDs
                # add them to a set since there may be duplicates, same rna-protien pair from different enrichments in different datasets
                for prot in proteinList:
                    pairsToWrite.add( transcriptID + "\t" + prot + "\n")


        for pair in pairsToWrite:                
            outFile2.write( pair)

        print "read_enrichment_data: read %s lines." % nLines
        print "read_enrichment_data: wrote %s transcript-protein pairs." % len( pairsToWrite)


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
                             help='Rainet database to be used.')
        parser.add_argument('enrichmentData', metavar='enrichmentData', type=str,
                             help='Data from enrichment analysis. This file should be pre-filtered for the wanted enrichments. Last column needs to be dataset of origin (e.g. WanCluster, BioplexCluster, etc), including change in header with "dataset" column.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
        parser.add_argument('--complexDatasets', metavar='complexDatasets', type=str, default = "WanCluster,BioplexCluster,CustomCluster",
                             help='List of datasets to withdrawn data from RAINET DB. Comma-separated.')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        instance = Enrichment2Protein( args.rainetDB, args.enrichmentData, args.outputFolder, args.complexDatasets)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        Timer.get_instance().step( "Read RainetDB table..")
        instance.read_rainet_db( )

        Timer.get_instance().step( "Read enrichment file..")            
        instance.read_enrichment_data()


        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

