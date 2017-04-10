
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
# 1) To separate well each complex, each protein in a complex is identified by a composed protein-complex key, so that the same protein in different complexes is interpreted differently
# 2) Using lncRNA gene IDs instead of transcript IDs
#===============================================================================

class PrepareDiseaseNetwork(object):


    NETWORK_OUTPUT_FILE = "/disease_network.tsv"
    NODE_INFO_OUTPUT_FILE = "/disease_network_node_info.tsv"
    CLUST_N_SEE_OUTPUT_FILE = "/disease_network_clustnsee.cns"

    def __init__(self, rainetDB, commonDiseaseFile, outputFolder):

        self.rainetDB = rainetDB
        self.commonDiseaseFile = commonDiseaseFile
        self.outputFolder = outputFolder

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
        # Output file, disease network
        #===============================================================================
        # example file
        # source node (protein or RNA) | target node (protein) | complexID (edge att) | (lncRNA) disease (edge att, if PRI)
        # prot1 | prot2 | complex1 | NA 
        # prot2 | prot1 | complex1 | NA
        # lncRNA1 | prot1 | complex1 | Telangiectasia
 
        # Note: only want information of complexes that are scaffolded by lncRNAs
        # Note: complexes will be represented as a clique, all proteins in a complex need to be interacting
        # Note: lncRNA-protein interactions come from the PRI table
        # Note: lncRNA used here is gene ID
         
        complexInfo = set( ) # key -> complex ID. Objective to only write info about a complex once
        sourceNodeDisease = set() # key -> proteinID, for proteins with disease
        sourceInteracting = set() # key -> proteinID, for RNA-interacting proteins
        doneInteractions = set() # prot1_prot2, to avoid duplication of interactions
        proteinNodes = set() # stores protein node IDs
        rnaNodes = set() # stores RNA node IDs
        proteinsPerComplex = {} # key -> complexID, value -> list of proteins
         
        networkInfo = [] # stores info that will be written in network output file
         
        # read file which contains one line per lncRNA-complex-protein_disease       
        inFile = open( self.commonDiseaseFile, "r")
        header = inFile.readline()
        for match in inFile:
            match = match.strip()
            spl = match.split("\t")
 
            #e.g.: transcriptID    GeneID  proteinID       wordsMatched    transcriptDisease       proteinDisease  complexID       complexProteins complexInteractions
            #e.g.: ['ENST00000539975', 'SNHG1', 'P00533', "'cell;lung;cancer'", 'non-small cell lung cancer', 'LUNG CANCER |  | ALVEOLAR CELL CARCINOMA, INCLUDED;; ADENOCARCINOMA OF LUNG, INCLUDED;; NONSMALL CELL LUNG CANCER, INCLUDED;; LUNG CANCER, PROTECTION AGAINST, INCLUDED', '313|NetworkModule', 'P31749,Q9Y4K3,P07900,Q13077,P00533,O15111,Q16543,A8K0R7,Q9P1W9,P21580,P04179,Q96F81,O75541,Q12933,P08238,P46778,Q92995,Q9H4A3,O14920,P63104,Q99558,Q9UN81,Q86YJ5,Q13571,P04049,P62258', 'P08238,P07900,O15111\n']
                         
            ##########################
            # Add transcript interaction info
            ##########################
            # lncRNA1 | prot1 | complex1 | Telangiectasia
            # Note: add a line for each lncRNA-protein_interacting in complex
 
            interactingProt = spl[8].strip().split(",")
 
            complexID = spl[6]
            geneID = spl[1]
            disease = spl[4]
             
            for protID in interactingProt:
                # composed key with protein and complex
                protTag = protID + "_" + complexID
                # add a line for each lncRNA-protein interaction
                networkInfo.append( [geneID, protTag, complexID, disease] )
                sourceInteracting.add( protTag)
            
            rnaNodes.add( geneID)
            
            diseaseProt = spl[2] + "_" + complexID
            sourceNodeDisease.add( diseaseProt)
             
            ##########################
            # Add protein complex interaction info (clique)
            ##########################           
            # prot1 | prot2 | complex1 | NA
 
            # only write complex interactions once for each complex            
            if complexID not in complexInfo:
                # get all proteins of complex
                complexProteins = spl[ 7].strip().split(",")
                 
                # for each pair of proteins, write interaction (one-way)
                for prot1 in complexProteins:
                    
                    # add proteins to complex structure
                    if complexID not in proteinsPerComplex:
                        proteinsPerComplex[ complexID] = []
                    proteinsPerComplex[ complexID].append( prot1 + "_" + complexID)
                    
                    for prot2 in complexProteins:
                        if prot1 != prot2: # avoid self interaction

                            # composite key
                            prot1Tag = prot1 + "_" + complexID
                            prot2Tag = prot2 + "_" + complexID
                            
                            if prot1 + "_" + prot2 not in doneInteractions: 
                                networkInfo.append( [prot1Tag, prot2Tag, complexID, "NA"] )                                
                                proteinNodes.add( prot1Tag)
                                proteinNodes.add( prot2Tag)
                            
                            # avoid repetition of inverse interaction
                            doneInteractions.add( prot1 + "_" + prot2)
                            doneInteractions.add( prot2 + "_" + prot1)
                         
            complexInfo.add( complexID)
 
        print len( proteinsPerComplex)
 
        print len( networkInfo)
 
        # Write network output file
        outFile = open( self.outputFolder + PrepareDiseaseNetwork.NETWORK_OUTPUT_FILE, "w")
        outFile.write( "Source_node\tTarget_node\tComplexID\tDisease_association\n")
        for line in networkInfo:
            outFile.write( "\t".join( line) + "\n")
                         
        outFile.close()
 
        print "read_common_disease_file: %s disease-associated complexes." %  len( complexInfo)
        # TODO: correct this, because does not work with composite keys
#        print "read_common_disease_file: %s disease-associated proteins." %  len( sourceNodeDisease)


        #===============================================================================
        # Output file, node info table
        #===============================================================================

        outFile2 = open(self.outputFolder + PrepareDiseaseNetwork.NODE_INFO_OUTPUT_FILE, "w")
        outFile2.write("Node\tType\tDisplay name\tDisease node\tRNA-interacting protein\n")

        # Write info of RNA nodes
        for rna in rnaNodes:
            outFile2.write("%s\tLncRNA\%s\t1\t1\n" % ( rna, rna) )

        for prot in proteinNodes:
            # tag disease proteins
            isDisease = 0
            if prot in sourceNodeDisease:
                isDisease = 1
                
            # tag interacting proteins
            isInteracting = 0
            if prot in sourceInteracting:
                isInteracting = 1

            # TODO: write protein name instead of uniprotAC
            outFile2.write("%s\tProtein\t%s\t%s\t%s\n" % ( prot, prot, isDisease, isInteracting ) )

        outFile2.close()


        #===============================================================================
        # Output file, clust N see partitions
        #===============================================================================
        # Note: each lncRNA in its own partition
        # Note: each complex in its own partition
        # Note: all nodes in partition needs also to be present in network

        # Example format      
        #         # ClustnSee analysis export
        #
        #         >ClusterID:1||
        #         prot1
        #
        #         >ClusterID:1||
        #         prot1
        #         prot2

        outFile3 = open( self.outputFolder + PrepareDiseaseNetwork.CLUST_N_SEE_OUTPUT_FILE, "w")
        outFile3.write( "#ClustnSee analysis export\n\n")

        # keep track of cluster ID, an incremental number
        clustID = 1

        # write lncRNA partitions
        for rna in rnaNodes:
            outFile3.write( ">ClusterID:%i||\n" % clustID)
            outFile3.write( rna + "\n")        
            outFile3.write("\n")
            clustID += 1

        # write complex partitions
        for comp in proteinsPerComplex:
            outFile3.write( ">ClusterID:%i||\n" % (clustID) )
            for prot in proteinsPerComplex[ comp]:
                outFile3.write( prot + "\n")                        
            outFile3.write("\n")
            clustID += 1

        outFile3.close()

        # TODO: unittest: ensure that at least one protein in each complex has disease


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
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        instance = PrepareDiseaseNetwork( args.rainetDB, args.commonDiseaseFile, args.outputFolder)
            
        instance.run()
            
        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())
