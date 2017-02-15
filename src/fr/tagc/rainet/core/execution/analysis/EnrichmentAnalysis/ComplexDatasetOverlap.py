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
from fr.tagc.rainet.core.data.InteractingProtein import InteractingProtein
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
# Started 15-Feb-2017
# Diogo Ribeiro
DESC_COMMENT = "Script to measure the overlap or redundancy between protein complex or module datasets."
SCRIPT_NAME = "ComplexDatasetOverlap.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read Rainet DB containing complex tables and interacting proteins
# 2) Calculate intra dataset overlap
# 3) Calculate inter dataset overlaps
# 4) Return text files for plotting in R
#===============================================================================

#===============================================================================
# Processing notes:
# 1)
# 2)
#===============================================================================

class ComplexDatasetOverlap(object):
    
    # correspondance between the base table and associated data
    ANNOTATION_TABLES_DICT = {"NetworkModule" : "ProteinNetworkModule",
                              "ReactomePathway" : "ProteinReactomeAnnotation",
                              "KEGGPathway" : "ProteinKEGGAnnotation",
                              "BioplexCluster" : "ProteinBioplexAnnotation",
                              "WanCluster" : "ProteinWanAnnotation",
                              "CorumCluster" : "ProteinCorumAnnotation",
                              "CustomCluster" : "ProteinCustomAnnotation"}
        
        
    DEFAULT_DATASET_LIST = "BioplexCluster,WanCluster,CorumCluster,CustomCluster,NetworkModule"

    OUTPUT_FILE_INTRA_DATASET = "intra_dataset_results.tsv"
    OUTPUT_FILE_INTER_DATASET = "inter_dataset_results.tsv"

    
    def __init__(self, rainetDB, outputFolder, useInteractingProteins, listDatasets):

        self.rainetDB = rainetDB
        self.outputFolder = outputFolder
        self.useInteractingProteins = useInteractingProteins

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)

        # Process list datasets
        self.listDatasets = []
        datasets = listDatasets.split(",")
        for dataset in datasets:
            if dataset in ComplexDatasetOverlap.ANNOTATION_TABLES_DICT:
                self.listDatasets.append( dataset)
            else:
                raise RainetException( "__init__: dataset not found: %s" % ( dataset ) )              
        

    # #
    # Get correspondence of complex to protein IDs from RainetDB, for all datasets
    def read_rainet_db(self):
        
        datasetDict = {} # key -> dataset name, value -> dict. key -> annotationID, value -> list of proteins in annotation
        
        for dataset in self.listDatasets:

            tableNameAnnotation = ComplexDatasetOverlap.ANNOTATION_TABLES_DICT[ dataset]
                 
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
     
            Logger.get_instance().info( "read_rainet_db: dataset: %s. read %s annotations" % ( dataset, len( complexAnnotDict ) ))
     
            datasetDict[ dataset] = complexAnnotDict

        
        Logger.get_instance().info( "read_rainet_db: processed %s datasets." % ( len( datasetDict)) )
    
        self.datasetDict = datasetDict
    
    # #
    # Calculates overlap of annotations within a dataset.
    def intra_dataset_overlap(self):

        #===================================================================   
        # Make all pairwise calculations
        #=================================================================== 
        
        intraDatasetOverlap = {} # key -> dataset, value -> dict. key -> annotationID, value -> list of overlaps
        
        for dataset in self.datasetDict:

            count = 0

            doneCalculations = set()
        
            if dataset not in intraDatasetOverlap:
                intraDatasetOverlap[ dataset] = {}
            
            datasetData = self.datasetDict[ dataset]
            
            for complex1 in datasetData:
                
                for complex2 in datasetData:
                
                    tag = complex1 + "|" + complex2
                    tag2 = complex2 + "|" + complex1
                    
                    # do not calculate overlap between itself
                    if complex1 == complex2:
                        continue
                    # skip in case calculation is already done, as it is done two-ways at once
                    if tag in doneCalculations or tag2 in doneCalculations:
                        continue

                    if complex1 not in intraDatasetOverlap[ dataset]:
                        intraDatasetOverlap[ dataset][ complex1] = []
      
                    if complex2 not in intraDatasetOverlap[ dataset]:
                        intraDatasetOverlap[ dataset][ complex2] = []
                    
                    perc1, perc2 = ComplexDatasetOverlap.group_overlap( datasetData[ complex1], datasetData[ complex2])

                    intraDatasetOverlap[ dataset][ complex1].append( perc1)
                    intraDatasetOverlap[ dataset][ complex2].append( perc2)

                    doneCalculations.add( tag)
                    doneCalculations.add( tag2)
                    count+=1


            # confirm that we are not calculating twice the same comparison
            assert( count * 2 == len( doneCalculations))

            Logger.get_instance().info( "intra_dataset_overlap: dataset: %s. %s overlaps calculated." % ( dataset, count ))

        #===================================================================   
        # Write output file
        #=================================================================== 
        
        outFile = open( self.outputFolder + ComplexDatasetOverlap.OUTPUT_FILE_INTRA_DATASET, "w")

        outFile.write( "dataset\tannotation_id\tmean_overlap\tall_annot_overlap\thigh_annot_overlap\n")
        
        intraDatasetOverlapResults = {}
        
        for dataset in intraDatasetOverlap:
            if dataset not in intraDatasetOverlapResults:
                intraDatasetOverlapResults[ dataset] = {}

            for annot in sorted( intraDatasetOverlap[ dataset]):
                
                if annot not in intraDatasetOverlapResults[ dataset]:
                    intraDatasetOverlapResults[ dataset][ annot] = ""
                
                # calculate the mean percentage of overlap for a complex, compared to all other complexes of the dataset
                meanOverlap = sum( intraDatasetOverlap[ dataset][ annot]) / len( intraDatasetOverlap[ dataset][ annot])

                # calculate in how many complexes there is any overlap
                anyOverlap = len( [x for x in intraDatasetOverlap[ dataset][ annot] if x != 0])

                # calculate in how many complexes there is overlap above 50%
                highOverlap = len( [x for x in intraDatasetOverlap[ dataset][ annot] if x > 50])

                outputText = "%s\t%s\t%.2f\t%s\t%s\n" % (dataset, annot, meanOverlap, anyOverlap, highOverlap)
                
                intraDatasetOverlapResults[ dataset][ annot] = outputText

                outFile.write( outputText)

        outFile.close()
        
        self.intraDatasetOverlap = intraDatasetOverlap
        self.intraDatasetOverlapResults = intraDatasetOverlapResults



    # #
    # Calculates percentage of group1 that is overlapped by group2, a vice versa
    @staticmethod
    def group_overlap( group1, group2):
        
        set1 = set( group1)
        set2 = set( group2)
        
        overlap = set1.intersection( set2)
        
        perc1 = len( overlap) * 100.0 / len( group1)
        perc2 = len( overlap) * 100.0 / len( group2)
        
        return perc1, perc2
                

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
                             help='Rainet database to be used. Needs to have protein annotation datasets and interactingProteins tables up-to-date.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Folder where output will be written.')
        # optional args
        parser.add_argument('--useInteractingProteins', metavar='useInteractingProteins', type=int, default = 1,
                             help='Whether to consider, for a complex, only proteins with interactions or all proteins in the complex. (Default = 1).')
        parser.add_argument('--listDatasets', metavar='listDatasets', type=int, default = "BioplexCluster,WanCluster,CorumCluster,CustomCluster,NetworkModule",
                             help='Which annotation datasets to use for this analysis. (Default = %s ).' % ComplexDatasetOverlap.DEFAULT_DATASET_LIST)
        

        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        currentRun = ComplexDatasetOverlap( args.rainetDB, args.outputFolder, args.useInteractingProteins, args.listDatasets)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        Timer.get_instance().step( "Read RainetDB table..")
        currentRun.read_rainet_db( )
        
        Timer.get_instance().step( "Calculate intra dataset overlap..")       
        currentRun.intra_dataset_overlap()

        #TODO: filter proteins by presence on interactingProteins table

        # Stop the chrono
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

