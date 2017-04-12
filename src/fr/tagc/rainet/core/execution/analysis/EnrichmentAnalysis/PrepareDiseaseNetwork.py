
import sys
import os
import argparse

# import numpy as np
# import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.RNA import RNA

#===============================================================================
# Started 10-April-2017
# Diogo Ribeiro
DESC_COMMENT = "Script to prepare files to produce a disease network."
SCRIPT_NAME = "PrepareDiseaseNetwork.py"
#===============================================================================


#===============================================================================
# General plan:
# 1) Read manually curated disease matching file, containing complexes and interacting proteins
# 2) Produce a file for producing a cytoscape network view and associated data tables
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Network at the lncRNA-gene and Complex level (not individual transcripts nor individual proteins)
#===============================================================================

class PrepareDiseaseNetwork(object):


    NETWORK_OUTPUT_FILE = "/disease_network.tsv"
    NODE_INFO_OUTPUT_FILE = "/disease_network_node_info.tsv"
#     CLUST_N_SEE_OUTPUT_FILE = "/disease_network_clustnsee.cns"

    def __init__(self, rainetDB, commonDiseaseFile, outputFolder, proteinShare):

        self.rainetDB = rainetDB
        self.commonDiseaseFile = commonDiseaseFile
        self.outputFolder = outputFolder
        self.proteinShare = proteinShare

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()
        

    # #
    # Function to read matching disease file and produce network file and associated tables
    def read_common_disease_file(self):
        
        #===============================================================================
        # Read input file
        #===============================================================================
        # Note: input file contains one line per lncRNA-complex-protein-disease
        # Note: a complex can have several diseases, but an RNA only one (or vice versa), therefore a lncRNA-complex enrichment should only contain the diseases in common between the two
        #
        # Example format    
        # transcriptID    GeneID  proteinID       wordsMatched    transcriptDisease       proteinDisease  complexID       complexProteins complexInteractions
        # ['ENST00000539975', 'SNHG1', 'P00533', "'cell;lung;cancer'", 'non-small cell lung cancer', 'LUNG CANCER |  | ALVEOLAR CELL CARCINOMA, INCLUDED;; ADENOCARCINOMA OF LUNG, INCLUDED;; NONSMALL CELL LUNG CANCER, INCLUDED;; LUNG CANCER, PROTECTION AGAINST, INCLUDED', '313|NetworkModule', 'P31749,Q9Y4K3,P07900,Q13077,P00533,O15111,Q16543,A8K0R7,Q9P1W9,P21580,P04179,Q96F81,O75541,Q12933,P08238,P46778,Q92995,Q9H4A3,O14920,P63104,Q99558,Q9UN81,Q86YJ5,Q13571,P04049,P62258', 'P08238,P07900,O15111\n']

        proteinNodes = set() # stores protein node IDs
        rnaNodes = set() # stores RNA node IDs
        proteinsPerComplex = {} # key -> complexID, value -> set of proteins
        diseasePerComplex = {} # key -> complexID, value -> set of diseases
        
        # stores info that will be written in network output file
        networkInfo = {} # key -> lncRNA + complex composite key, value -> info for that association
        networkInfoDisease = {} # key -> lncRNA + complex composite key, value -> diseases for that association
        
        inFile = open( self.commonDiseaseFile, "r")
        header = inFile.readline()
        for match in inFile:
            match = match.strip()
            spl = match.split("\t")
                          
            ##########################
            # Add transcript interaction info
            ##########################
            # lnc1 | complex1 | 5 | pri 
            # Note: add a line for each lncRNA-complex interaction
            # Note: if there is several diseases being shared, these will be acumulated in the 'disease' column, not in another entry
 
            complexID = spl[6]
            geneID = spl[1]
            disease = spl[4]
            interactingProt = spl[8].strip().split(",")
            
            compositeKey = geneID + "_" + complexID

            # add basic interaction info  
            if compositeKey not in networkInfo:
                networkInfo[ compositeKey] = []
            networkInfo[ compositeKey] = [ geneID, complexID, str( len( interactingProt)), "pri"]

            # add disease association info to edges
            if compositeKey not in networkInfoDisease:
                networkInfoDisease[ compositeKey] = set()
            networkInfoDisease[ compositeKey].add( disease)

            # add proteins to complex structure
            if complexID not in proteinsPerComplex:
                proteinsPerComplex[ complexID] = set()
                complexProteins = spl[ 7].strip().split(",")

                for protein in complexProteins:
                    proteinsPerComplex[ complexID].add( protein)
                    proteinNodes.add( protein)

            # add disease information to complexes
            if complexID not in diseasePerComplex:
                diseasePerComplex[ complexID] = set()
            diseasePerComplex[ complexID].add( disease)
            
            rnaNodes.add( geneID)

        assert len( networkInfo) == len( networkInfoDisease)
             
        print "read_common_disease_file : number of complex nodes: %s " % len( proteinsPerComplex)
        print "read_common_disease_file : total number of proteins: %s " % len( proteinNodes)
        print "read_common_disease_file : number of RNA nodes: %s " % len( rnaNodes)
        print "read_common_disease_file : number of RNA-complex interactions (edges): %s " % len( networkInfo)



        #===============================================================================
        # Output file, disease network
        #===============================================================================
        # example file
        # source node (complex or RNA) | target node (complex) | number partners | interaction type | diseases
        # lnc1 | complex1 | 5 | pri | breast cancer, skin cancer
        # complex1 | complex2 | 2 | protein share | NA
 
        # Note: only want information of complexes that are scaffolded by lncRNAs
        # Note: lncRNA used here is gene ID
 
        # Write PRI
        outFile = open( self.outputFolder + PrepareDiseaseNetwork.NETWORK_OUTPUT_FILE, "w")
        outFile.write( "Source_node\tTarget_node\tNumber partners\tInteraction type\tDisease associations\n")
        for assoc in sorted( networkInfo):
            # add disease info to network line
            diseaseMerge = ", ".join( networkInfoDisease[ assoc])
            networkInfo[ assoc].append( diseaseMerge )
            outFile.write( "\t".join( networkInfo[ assoc]) + "\n")

        # Write complex protein share
        # complex1 | complex2 | 2 | protein share | NA

        complexShare = {}
        if self.proteinShare != -1:
            complexShare = self.calculate_protein_share(proteinsPerComplex)
            
            for pair in complexShare:
                comp1, comp2 = pair.split( "_")
                outFile.write( "%s\t%s\t%i\tprotein share\t'\n" % ( comp1, comp2, complexShare[ pair] ) )
                    
        outFile.close()
        
        #===============================================================================
        # Output file, node info table
        #===============================================================================
        # Node ID | Node type | Complex size (if applicable) | Complex Disease
 
        outFile2 = open(self.outputFolder + PrepareDiseaseNetwork.NODE_INFO_OUTPUT_FILE, "w")
        outFile2.write("Node\tType\tComplex size\tComplex disease\n")
 
        # Write info of RNA nodes
        for rna in rnaNodes:
            # 50 is default size for an RNA node
            outFile2.write( "%s\tLncRNA\t50\t%s\n" % ( rna, rna) )
 
        # Write info for complex nodes
        for comp in proteinsPerComplex:
            outFile2.write( "%s\tComplex\t%i\t%s\n" % ( comp, len( proteinsPerComplex[ comp]), ", ".join( diseasePerComplex[ comp]) ) )
 
        outFile2.close()

        self.networkInfo = networkInfo
        self.complexShare = complexShare


    # #
    # Calculate number of proteins shared between all combinations of two complexes (one-way, non-self)
    def calculate_protein_share(self, proteinsPerComplex):
        
        #===============================================================================
        # Calculate complex protein share
        #===============================================================================
        # Note: complex interactions are number of proteins being shared between them
        # Note: only write interactions where there is at least one protein being shared

        complexShare = {} # key -> e.g. complex1_complex2, value -> number proteins being shared
    
        donePairs = set() # store done comparisons, so that we only do them one way
        for comp1 in proteinsPerComplex:
            for comp2 in proteinsPerComplex:
                if comp1 != comp2:
                    
                    tag1 = comp1 + "_" + comp2
                    tag2 = comp2 + "_" + comp1

                    if tag1 not in donePairs:
                        overlap = len( proteinsPerComplex[ comp1].intersection( proteinsPerComplex[ comp2]) )

                        if overlap >= self.proteinShare:
                            complexShare[ tag1] = overlap
                    
                    donePairs.add( tag1)
                    donePairs.add( tag2)
        
        possibleCombinations = ((len( proteinsPerComplex)-1) * len( proteinsPerComplex)) / 2
        assert( len( complexShare) <= possibleCombinations)
        
        print "calculate_protein_share : number of complex combinations sharing proteins: %s, out of %s possible combinations" % ( len( complexShare), possibleCombinations )
        
        return complexShare


    def run(self):
        
        self.read_common_disease_file()



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
        parser.add_argument('commonDiseaseFile', metavar='commonDiseaseFile', type=str,
                             help='File with matches between lncRNA-complex diseases, the manually curated output from CommonLncRNAProteinDisease.py.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
        parser.add_argument('--proteinShare', metavar='--proteinShare', type=int, default = -1,
                             help='Whether to add interactions for protein share between complexes, and how many minimum proteins needed for writing interaction. Default = -1 (OFF)')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        instance = PrepareDiseaseNetwork( args.rainetDB, args.commonDiseaseFile, args.outputFolder, args.proteinShare)           
            
        instance.run()
        
        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())
