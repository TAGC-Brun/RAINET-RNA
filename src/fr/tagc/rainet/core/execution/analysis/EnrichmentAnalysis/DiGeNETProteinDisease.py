
import argparse
import os
import sys

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
from fr.tagc.rainet.core.util.time.Timer import Timer

# from fr.tagc.rainet.core.data.Protein import Protein
# from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference


# import numpy as np
# import pandas as pd
# from fr.tagc.rainet.core.util.data.DataManager import DataManager
#===============================================================================
# Started 10-January-2017
# Diogo Ribeiro
DESC_COMMENT = "Script to process and reformat DiGeNET."
SCRIPT_NAME = "DiGeNETProteinDisease.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read DiGeNET data, map protein IDs
# 2) Return file with protein-disease correspondence
#===============================================================================

#===============================================================================
# Processing notes:
#===============================================================================

class DiGeNETProteinDisease(object):

    # Constants
        
    OUTPUT_FILE = "/DiGeNET_protein_disease_description.txt"
        
    def __init__(self, digenetIDMappingFile, digenetFile, outputFolder):

        self.digenetIDMappingFile = digenetIDMappingFile
        self.digenetFile = digenetFile
        self.outputFolder = outputFolder

#         # Build a SQL session to DB
#         SQLManager.get_instance().set_DBpath(self.rainetDB)
#         self.sql_session = SQLManager.get_instance().get_session()


        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)
            
            
#     # #
#     # Get correspondence of protein IDs from RainetDB
#     def read_rainet_db(self):
# 
#         proteinCrossReference = {} # key -> uniprotID, val -> uniprotAC
# 
#         #===============================================================================
#         # Query the protein table
#         #===============================================================================
# 
#         query = self.sql_session.query( Protein.uniprotAC, Protein.uniprotID ).all()
# 
#         # uniprotAC example: "P31946"
#         # uniprotID example: "1433B_HUMAN"
# 
#         setOfProteinIDs = set()
#         for uniprotAC, uniprotID in query:
#             
#             # remove the "_HUMAN" tag to the ID
#             if "_HUMAN" in uniprotID:
#                 uniprotID = uniprotID.split("_HUMAN")[0]
#             
#             if uniprotID not in proteinCrossReference:
#                 proteinCrossReference[ uniprotID] = set( )
#                 
#             proteinCrossReference[ uniprotID].add( str( uniprotAC) )
# 
#             setOfProteinIDs.add( uniprotAC)
# 
#         print "read_rainet_db: Read Protein table."
#         print "read %s entries" % len( query)
#         print "read %s IDs mapped for %s proteins" % (len( proteinCrossReference), len(setOfProteinIDs) )
# 
#         #===============================================================================
#         # Query the proteinCrossReference table
#         #===============================================================================
# 
#         query = self.sql_session.query( ProteinCrossReference.protein_id, ProteinCrossReference.crossReferenceID ).all()
# 
#         # query ProteinCrossReferences on GeneCards database only
# #        query = self.sql_session.query( ProteinCrossReference.protein_id, ProteinCrossReference.crossReferenceID ).filter(ProteinCrossReference.sourceDB=="GeneCards").all()
# 
#         for uniprotAC, crossID in query:
#             
#             if crossID not in proteinCrossReference:
#                 proteinCrossReference[ crossID] = set( )
#                 
#             proteinCrossReference[ crossID].add( str( uniprotAC) )
# 
#             setOfProteinIDs.add( uniprotAC)
# 
#         print "read_rainet_db: Read ProteinCrossReference table."
#         print "read %s entries" % len( query)
#         print "read %s IDs mapped for %s proteins" % (len( proteinCrossReference), len(setOfProteinIDs) )
# 
#         self.proteinCrossReference = proteinCrossReference

    # #
    # Read ID mapping file provided by DiGeNET mapa_geneid_4_uniprot_crossref.tsv
    def read_digenet_id_mapping_file( self): 
        
        #===============================================================================
        # Read DiGeNET mapa_geneid_4_uniprot_crossref.tsv
        #===============================================================================
        # Example format:
        # UniProtKB       GENEID
        # P04217  1
        # P01023  2

        entrezUniprotDict = {} # key -> entrez ID, value -> uniprotAC
        setOfProteinIDs = set()
        nLines = 0

        with open( self.digenetIDMappingFile, "r") as inFile:
            header = inFile.readline()
            nLines += 1
            
            for line in inFile:
                spl = line.strip().split( "\t")

                nLines += 1
                
                uniprotAC = spl[0]
                entrez = spl[1]
        
                if entrez not in entrezUniprotDict:
                    entrezUniprotDict[ entrez] = set()
                    
                entrezUniprotDict[ entrez].add( uniprotAC)

                setOfProteinIDs.add( uniprotAC)

        print "read_digenet_id_mapping_file: Read DiGeNET ID mapping file."
        print "read_digenet_id_mapping_file: read %s entries" % nLines
        print "read_digenet_id_mapping_file: read %s IDs mapped for %s proteins" % (len( entrezUniprotDict), len(setOfProteinIDs) )

        self.entrezUniprotDict = entrezUniprotDict


    # #
    # Read curated_gene_disease_associations.txt file from DiGeNET database and map IDs, write output file
    def read_curated_digenet_file( self):

        #===============================================================================
        # Read DiGeNET curated_gene_disease_associations.txt
        #===============================================================================    
        # Example format:    
        # diseaseId       geneId  score   geneName        description     diseaseName     sourceId        NofPmids        NofSnps
        # umls:C0000737   3440    0.12    IFNA2   interferon, alpha 2     Abdominal Pain  CTD_human       1       0
        # umls:C0000744   4547    0.490867606404876       MTTP    microsomal triglyceride transfer protein        Abetalipoproteinemia    CLINVAR,CTD_human,ORPHANET,UNIPROT      33      5
        
        # Note: commented lines start with "#"


        #===============================================================================
        # Output file
        #===============================================================================
        # proteinID\tdisease_ID\tdisease_name

        outFile = open( self.outputFolder + DiGeNETProteinDisease.OUTPUT_FILE, "w")

        notFound = set()
        found = set()
        nLines = 0

        with open( self.digenetFile, "r") as inFile:
            for line in inFile:

                nLines += 1

                # skip commented lines
                if line.startswith( "#"):
                    continue

                if line.startswith( "diseaseId"):
                    header = line
                    continue

                spl = line.strip().split( "\t")

                # retrieve wanted attributes                
                diseaseID = spl[0]
                geneID = spl[1]
                score = spl[2]
                diseaseName = spl[5]
                
                if geneID in self.entrezUniprotDict:
                    found.add( geneID)
                    # retrieve all uniprotACs associated to this gene name / protein ID
                    uniprotAC = self.entrezUniprotDict[ geneID]
                    # write output line for each matched uniprotAC
                    for uni in uniprotAC:
                        outFile.write( "%s\t%s\t%s\n" % ( uni, diseaseID, diseaseName) )
                else:
                    notFound.add( geneID)

        print "read_curated_digenet_file: read %s lines." % nLines
        print "read_curated_digenet_file: could not map %s DiGeNET protein IDs to uniprotAC." % len( notFound)
        print "read_curated_digenet_file: mapped %s DiGeNET protein IDs to uniprotAC." % len( found)
 
        # Note: curated file seems to have a minimum of 0.12 as score, maximum of 0.86.

    # #
    # Read all_gene_disease_associations.tsv file from DiGeNET database and map IDs, write output file
    def read_all_gene_disease_digenet_file( self):
        
        #TODO: ..
        pass


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
#         parser.add_argument('rainetDB', metavar='rainetDB', type=str,
#                              help='Rainet database to be used.')
        parser.add_argument('digenetIDMappingFile', metavar='digenetIDMappingFile', type=str,
                             help='Path to DiGeNET file to be used.')
        parser.add_argument('digenetFile', metavar='digenetFile', type=str,
                             help='Path to DiGeNET file to be used.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        instance = DiGeNETProteinDisease( args.digenetIDMappingFile, args.digenetFile, args.outputFolder)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        Timer.get_instance().step( "Read digenet ID Mapping File ..")
        instance.read_digenet_id_mapping_file( )

        Timer.get_instance().step( "Read DiGeNET file..")            
        instance.read_curated_digenet_file()
        
        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

