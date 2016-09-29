import sys
import os
import argparse
import igraph

import numpy as np
# import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

# from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.data.DataManager import DataManager

from fr.tagc.rainet.core.data.Protein import Protein

#===============================================================================
# Started 28-Sep-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to analyse properties of RNA-protein interactions on the protein-protein interaction network."
SCRIPT_NAME = "NetworkScoreAnalysis.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read PPI network, calculate degree of each protein
# 2) Read catRAPID interactions file, pick top X protein binders, calculate several metrics of this set of proteins on the PPI
# 3) Compare real results against random control
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Using uniprot ID instead of uniprot AC
# 2) Ignoring interactions with proteins not in the PPI network. the fact that a protein is not in the PPI network, does not mean that protein has no known interactions, but that we can't easily apply any metrics, therefore we ignore those cases
#===============================================================================

class NetworkScoreAnalysis(object):
    
    #===============================================================================
    # NetworkScoreAnalysis Constants
    #===============================================================================        

    #===================================================================
    # Report files constants       
    #===================================================================

    PARAMETERS_LOG = "parameters.log"

#     # Annotation report
#     REPORT_PROT_PER_ANNOTATION = "prot_per_annotation.tsv"

       
    def __init__(self, networkFile, catrapidFile, rainetDBFile, topPartners, outputFolder, numberRandomizations):

        self.networkFile = networkFile
        self.catrapidFile = catrapidFile
        self.rainetDBFile = rainetDBFile
        self.topPartners = topPartners
        self.outputFolder = outputFolder
        self.numberRandomizations = numberRandomizations

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDBFile)
        self.sql_session = SQLManager.get_instance().get_session()

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)


#     # #
#     # Use RAINET DB to retrieve Protein cross references
#     def protein_cross_references(self):
# 
# 
#         proteinIDMapping = "proteinIDMapping"
# 
#         DataManager.get_instance().perform_query( proteinIDMapping, "query( Protein.uniprotAC, Protein.uniprotID).all()") 
# 
#         # Convert query into a dictionary
#         DataManager.get_instance().query_to_dict( proteinIDMapping, 1, 0)
#         proteinIDMappingDict = DataManager.get_instance().get_data( proteinIDMapping) # key -> uniprotID (e.g. .._HUMAN), val -> uniprotAC (e.g. P35670)
#         # formatting
#         proteinIDMappingDict = { prot : proteinIDMappingDict[prot][0] for prot in proteinIDMappingDict}
# 
#         self.proteinIDMappingDict = proteinIDMappingDict
   

    # #
    # Function to read binary PPI network file into an igraph object
    def read_network_file(self):

        listOfTuples = [] #e.g. [ (0,1), (1,2), (0,2) ]
        listOfNames = [] #e.g. [('PRRT3_HUMAN', 'TMM17_HUMAN')]
        
        dictNames = {} # key -> protein name, value -> graph index (internal id)
        idCounter = 0
        
        with open( self.networkFile, "r") as inFile:
            for line in inFile:
                line = line.strip()
                node1, node2 = line.split("\t")
        
                if node1 not in dictNames:
                    idCounter += 1
                    dictNames[ node1] = idCounter
        
                if node2 not in dictNames:
                    idCounter += 1
                    dictNames[ node2] = idCounter
        
                names = (node1, node2)
                ids = (dictNames[ node1], dictNames[ node2])
                
                listOfNames.append( names)
                listOfTuples.append( ids)

        assert( len( listOfNames) == len( listOfTuples) )
        
        graph = igraph.Graph(listOfTuples)
        
        Logger.get_instance().info( "NetworkScoreAnalysis.read_network_file : Number of PPIs %s" % ( len( listOfNames)) )
        Logger.get_instance().info( "NetworkScoreAnalysis.read_network_file : Number of unique proteins %s" % ( len( dictNames)) )
    
        # swap dictionary to have graph index as keys (note that there is 1-to-1 correspondence)
        proteinGraphIDDict = { dictNames[ prot] : prot for prot in dictNames} # key -> graph index (internal id), value -> protein name
        
        assert( len( dictNames) == len( dictNames))
    
        self.graph = graph
        self.proteinGraphIDDict = proteinGraphIDDict
        self.dictNames = dictNames
                
        # for testing purposes
        return graph, listOfNames, listOfTuples, dictNames
    
    
    # #
    # Function to calculate degree level for each protein in PPI network
    def calculate_protein_degree(self):

        degreeDict = {} # key -> uniprotID, val -> degree of protein
        
        g = self.graph

        proteinNotFound = set()

        # for each vertice in graph
        for v in g.vs:
            
            # look for protein name using index
            try:
                proteinName = self.proteinGraphIDDict[ v.index]
            except KeyError:
                if v.index == 0:
                    continue
                else:
                    raise RainetException( "NetworkScoreAnalysis.calculate_protein_degree : Could not find protein name for graph index %s" % ( v.index ) )
                    
#             try:
#                 proteinUniprotAC = self.proteinIDMappingDict[ proteinName]
#             except KeyError:
#                 Logger.get_instance().warning( "NetworkScoreAnalysis.calculate_protein_degree : Could not find protein uniprotac for protein name %s. We will be converting here to uniprotAC" % ( proteinName ) )
#                 
#                 if "_HUMAN" in proteinName:
#                     proteinUniprotAC = proteinName.split("_")[0]
#                 else:
#                     proteinNotFound.add( proteinName)
#                     continue
                        
            if proteinName not in degreeDict:
                degreeDict[ proteinName] = v.degree()
            else:
                raise RainetException( "NetworkScoreAnalysis.calculate_protein_degree : duplicate uniprotAC %s" % ( proteinName ) )
 
        Logger.get_instance().info( "NetworkScoreAnalysis.calculate_protein_degree : Could not find proteins uniprotac for protein names %s. These were discarded." % ( proteinNotFound ) )
 
        self.degreeDict = degreeDict


    # #
    # Read catrapid file, build dictionary for each RNA containing scores and interacting proteins.
    # Use uniprotID (e.g. .._HUMAN)
    def read_catrapid_file( self):

        # From template of ReadCatrapid.py

        #=======================================================================
        # Example file
        # sp|Q96DC8|ECHD3_HUMAN ENST00000579524   -12.33  0.10    0.00
        # sp|P10645|CMGA_HUMAN ENST00000516610    10.66   0.32    0.00
        # protein and rna separated by " ", other values separated by "\t"
        # 
        # Protein is always on left side, RNA in the right side.
        # Assumption that there only one interaction between each Protein-RNA pair
        #=======================================================================

        rnaTargets = {} # key -> transcript ID (ensembl..), value -> score

        allProtSet = set()
        allRNASet = set()
    
        #=======================================================================
        # read file
        #=======================================================================
        with open( self.catrapidFile, "r") as inFile:
            for line in inFile:

                spl = line.split(" ")
                
                protID = spl[0].split( "|")[2] #pick uniprotID instead of uniprotAC
                spl2 = spl[1].split( "\t")
                rnaID = spl2[0]
                score = float( spl2[1])
                
                # decided not to round score values, but may be useful if input file is too large
                scoreRounded = score
                #scoreRounded = round( score, 1) 
                           
                allRNASet.add( rnaID)
                allProtSet.add( protID)

                ## RNA side
                if rnaID not in rnaTargets:
                    rnaTargets[ rnaID] = {}

                if scoreRounded not in rnaTargets[ rnaID]:
                    rnaTargets[ rnaID][ scoreRounded] = set()
                rnaTargets[ rnaID][ scoreRounded].add( protID)

        
        assert( len(allRNASet) == len( rnaTargets))

        Logger.get_instance().info( "NetworkScoreAnalysis.read_catrapid_file : Number of unique RNAs = %s" % ( len(allRNASet) ) )
        Logger.get_instance().info( "NetworkScoreAnalysis.read_catrapid_file : Number of unique Proteins = %s" % ( len(allProtSet) ) )

        self.rnaTargets = rnaTargets
        self.allProtSet = allProtSet
        self.allRNASet = allRNASet


    # #
    # Pick top X proteins for each RNA, based on best score, in case of same score, pick first appearence in input file.
    def pick_top_proteins(self):

        rnaTargets = self.rnaTargets
        
        #=======================================================================
        # pick top X proteins for each RNA
        #=======================================================================
        
        rnaTops = {} # key -> transcript ID, value -> list of Top X proteins
        
        for rna in rnaTargets:
                        
            if rna not in rnaTops:
                rnaTops[ rna] = []              
            
            # switch of whether to keep searching for more scores to get more top proteins
            boo = 1
            
            sortedScores = sorted(rnaTargets[ rna].keys(), reverse = True)
            
            for score in sortedScores:
                if boo:
                    for prot in rnaTargets[ rna][ score]:

                        # add top proteins to rnaTops if top limit is not yet reached
                        if len( rnaTops[ rna]) < self.topPartners:
                            # check if protein present in provided graph, otherwise pick next protein
                            if prot in self.dictNames:
                                rnaTops[ rna].append( prot)
                            else:
                                # the fact that a protein is not in the PPI network, does not mean that protein has no known interactions, but that we can't easily apply any metrics, therefore we ignore those cases
                                Logger.get_instance().warning( "NetworkScoreAnalysis.pick_top_proteins : %s. 'Top' protein %s is not present in PPI network, will use next top protein." % ( rna, prot ) )
                                
                        else:
                            # if top is full, stop searching for more proteins
                            boo = 0
                            continue

            # Warn if there is not enough proteins to fill top
            if len( rnaTops[ rna]) < self.topPartners:
                Logger.get_instance().warning( "NetworkScoreAnalysis.pick_top_proteins : %s does not have enough interactions to fill provided top. %s proteins are used." % ( rna, len( rnaTops[ rna]) ) )

        self.rnaTops = rnaTops

        
    
    # #
    # For each RNA, calculate several metrics for their top protein partners in their PPI network
    def calculate_metrics(self):
                
        dictNames = self.dictNames    
        
        rnaTops = self.rnaTops        
        graph = self.graph
 
        #=======================================================================
        # Calculate mean of mean shortest path between top proteins
        #=======================================================================
        # for each RNA, for each top protein, calculate shortest path against each other top protein. Calculate mean for each protein, and then mean for each RNA.
 
        rnaShortestPath = {} # key -> transcript ID, val -> mean of mean shortest paths
 
        for rna in rnaTops:

            allIdx = [ dictNames[ prot] for prot in rnaTops[ rna]]
            
            assert( len( allIdx) == self.topPartners)

            meanShortestPaths = []

            for idx in allIdx:
                # get list of indexes other than current
                otherIdx = [i for i in allIdx if i != idx ] 
                # calculate shortest path between current node against all others
                shortestPaths = graph.shortest_paths( idx, otherIdx, mode = "OUT")[0]
                # calculate mean shortest path for a node
                meanShortestPath = np.mean( shortestPaths)
                meanShortestPaths.append( meanShortestPath)
            
            # calculate mean shortest path for current RNA
            meanRNAShortestPath = np.mean( meanShortestPaths)
            
            rnaShortestPath[ rna] = "%.2f" % meanRNAShortestPath

        
        print rnaShortestPath
         
             
        # TODO: when making random control, use same amount of proteins as existing for each RNA
 
     
     
        # takes time to run...
        #print (graph.average_path_length())

    
    
    # #
    # Central function to run functions in order
    def run( self):

        Logger.get_instance().info( "NetworkScoreAnalysis.run: Starting..." )

        #===================================================================
        # Initialising datasets
        #===================================================================

        Timer.get_instance().step( "Reading network file.." )        
        
        #self.protein_cross_references()
        self.read_network_file()
        self.calculate_protein_degree()

        Timer.get_instance().step( "Reading catrapid file.." )        

        self.read_catrapid_file()        
        self.pick_top_proteins()

         
#         Logger.get_instance().info( "NetworkScoreAnalysis.run : ... %s" % ( len( self.testContainer)) )
   

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
        parser.add_argument('networkFile', metavar='networkFile', type=str,
                             help='File with binary protein-protein interaction network. E.g. PRRT3_HUMAN\tTMM17_HUMAN.')
        parser.add_argument('catrapidFile', metavar='catrapidFile', type=str,
                             help='File with all vs all catrapid results. All interactions there will be processed. Ideally RNA and protein set would be filtered, but not at the cutoff/expression level.')
        parser.add_argument('rainetDBFile', metavar='rainetDBFile', type=str,
                             help='File with a RAINET Database to use for protein ID mapping.')
        parser.add_argument('topPartners', metavar='topPartners', type=int,
                             help='Number of top protein partners to consider for analysis, for each transcript.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=int,
                             help='Where to write output files.')
        # optional args
        parser.add_argument('--numberRandomizations', metavar='numberRandomizations', type=int, default = 1000,
                             help='Number of randomizations to be performed to calculate empirical p-value for each metric, for each transcript.')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        networkScoreAnalysis = NetworkScoreAnalysis( args.networkFile, args.catrapidFile, args.rainetDBFile, args.topPartners, args.outputFolder, args.numberRandomizations)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        networkScoreAnalysis.run( )

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

