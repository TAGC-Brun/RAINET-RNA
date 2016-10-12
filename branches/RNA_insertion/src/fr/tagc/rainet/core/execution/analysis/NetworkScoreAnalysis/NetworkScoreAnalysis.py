import sys
import os
import argparse
import igraph
import copy
import random

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
#
# Lionel/Christine neighbour score metric (LCneighbours):
# 
# * Let's consider the first N protein in the top list of interaction to a given RNA.
# 
# * For each of those proteins, locate the corresponding node in the network
# 
# * Compute the shortest distances from this node to the nodes of the (N-1) other proteins
# 
# * For each distance, if it is 1, add 1 to the score of the initial node, if it's 2, add 1/2 and if it's three add 1/3. Add 0 if the distance in greater than 3.
# 
# * Sum up the scores of the N proteins to get a global score and write it down
# 
# * Finally, compute this score varying N from 3 to 20 for instance...
# 
# * Display the evolution of the score on a graph....
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Using uniprot ID instead of uniprot AC
# 2) Ignoring interactions with proteins not in the PPI network. the fact that a protein is not in the PPI network, does not mean that protein has no known interactions, but that we can't easily apply any metrics, therefore we ignore those cases
# 3) For randomization: if there are no more proteins to sample with same degree, take proteins from closest degree possible.
# 4) If transcript has less interactions than number of wanted top proteins, transcript is not processed.
# 5) Potential issue: for picking random proteins we use whole PPI network, regardless of having interactions or not. This may be different protein set as we only have interactions for proteins <750aa
# 6) Assortativity: using undirected assortativity of graph, populating graph with presence/absence (1/0) of interaction based on a cutoff. If no interaction information (e.g. protein >750aa), this protein is treated as having no interaction (0).
#===============================================================================


class NetworkScoreAnalysis(object):
    
    #===============================================================================
    # NetworkScoreAnalysis Constants
    #===============================================================================        

    # Lionel metric scores, if shortest path is one, it will have a specific score, if shortest path value is not in the dictionary the "DISTANCE_ABOVE_LIMIT_SCORE" is used
    DISTANCE_SCORE = {
                      1 : 1,
                      2 : 0.5,
                      3 : 0.25 }

    DISTANCE_ABOVE_LIMIT_SCORE = -0.1

    # Number of interactions above this value will use too much memory
    MAXIMUM_NUMBER_VIABLE_INTERACTIONS = 30000000    

    # For assortativity metric: interaction propensity cutoff to distinguish interacting from non-interacting.
    SCORE_CUTOFF = 15

    #===================================================================
    # Report files constants       
    #===================================================================

    PARAMETERS_LOG = "parameters.log"

    # Metrics report
    REPORT_METRICS_OUTPUT = "metrics_per_rna.tsv"

       
    def __init__(self, networkFile, catrapidFile, topPartners, outputFolder, numberRandomizations):

        self.networkFile = networkFile
        self.catrapidFile = catrapidFile
        self.topPartners = topPartners
        self.outputFolder = outputFolder
        self.numberRandomizations = numberRandomizations

#         # Build a SQL session to DB
#         SQLManager.get_instance().set_DBpath(self.rainetDBFile)
#         self.sql_session = SQLManager.get_instance().get_session()

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)


        self.calculatedShortestPaths = {}


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
        
        proteinsPerDegreeDict = {} # key -> degree, val -> set of proteins with that degree
        
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
 
        # create dictionary with list of proteins per degree
        for prot in degreeDict:
            degree = degreeDict[ prot]
            if degree not in proteinsPerDegreeDict:
                proteinsPerDegreeDict[ degree] = set()
            proteinsPerDegreeDict[ degree].add( prot)
            

        ## list of degrees available, so that closest degree to the degree of interest can be probed
        #sortedDegrees = sorted( proteinsPerDegreeDict.keys()) 

        self.degreeDict = degreeDict
        self.proteinsPerDegreeDict = proteinsPerDegreeDict
        #self.sortedDegrees = sortedDegrees


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
        
        nlines = 0
        
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

                nlines += 1

                if nlines > NetworkScoreAnalysis.MAXIMUM_NUMBER_VIABLE_INTERACTIONS:
                    raise RainetException( "NetworkScoreAnalysis.read_catrapid_file : number of interactions is too large to be computable: %s interactions" % nlines)

        
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
            
            # store proteins that are in top interactors but not present in network, and therefore skipped in analysis
            skippedProteins = set()     
                   
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
                                skippedProteins.add( prot)                                
                        else:
                            # if top is full, stop searching for more proteins
                            boo = 0
                            continue

            # Warn about how many top proteins skipped because they are not present in PPI network
            # the fact that a protein is not in the PPI network, does not mean that protein has no known interactions, but that we can't easily apply any metrics, therefore we ignore those cases
            Logger.get_instance().warning( "NetworkScoreAnalysis.pick_top_proteins : %s. 'Top' proteins skipped for not being present in PPI network: %s. List: %s" % ( rna, len( skippedProteins), ",".join( skippedProteins) ) )

            # Warn if there is not enough proteins to fill top
            # transcript is not analysed as it wouldn't be comparable to others
            if len( rnaTops[ rna]) < self.topPartners:
                Logger.get_instance().warning( "NetworkScoreAnalysis.pick_top_proteins : %s does not have enough interactions (has %s) to fill provided top. No results for this transcript." % ( rna, len( rnaTops[ rna]) ) )
                del rnaTops[ rna]

        self.rnaTops = rnaTops

        
    # #
    # For each RNA, calculate several metrics for their top protein partners in their PPI network
    def calculate_metrics(self):
                        
        rnaTops = self.rnaTops        
 
        #=======================================================================
        # - Calculate mean of mean shortest path between top proteins
        # - Calculate Lionel metric (see details at start of script)
        # - Calculate assortativity
        #=======================================================================
        # Shortest path mean: for each RNA, for each top protein, calculate shortest path against each other top protein. Calculate mean for each protein, and then mean for each RNA.
        # LCneighbour metric: see top of script
    
        #=======================================================================
        # Write output file
        #=======================================================================

        outFile = open( self.outputFolder + "/" + NetworkScoreAnalysis.REPORT_METRICS_OUTPUT, "w" )

        # write header
        outFile.write("transcriptID\tLCneighbours\tLCneighboursRandom\tLCneighboursPval\tShortestPath\tShortestPathRandom\tShortestPathPval\tAssortativity\tAssortativityRandom\tAssortativityPval\n")

        rnaShortestPath = {} # key -> transcript ID, val -> mean of mean shortest paths
        rnaShortestPathRandom = {} # key -> transcript ID, val -> list of mean of mean shortest paths, one for each randomization
        rnaShortestPathPval = {} # key -> transcript ID, val -> pval
 
        lionelMetrics = {} # key -> transcript ID, val -> lionel metric
        lionelMetricsRandom = {} # key -> transcript ID, val -> list of lionel metrics, one for each randomization
        lionelMetricsPval = {} # key -> transcript ID, val -> pval
 
        assortativityMetric = {} # key -> transcript ID, val -> assortativity
        assortativityMetricRandom = {} # key -> transcript ID, val -> list of assortativity values, one for each randomization
        assortativityMetricPval = {} # key -> transcript ID, val -> pval
 
        count = 0
 
        ### calculate metrics for each RNA
        for rna in rnaTops:
                  
            count+=1
            if count % 100 == 0:
                Logger.get_instance().info( "NetworkScoreAnalysis.calculate_metrics: processed %s transcripts.." % ( count) )           
        
            ## calculate metrics for real data

            topProteins = rnaTops[ rna]

            meanRNAShortestPath, lionelMetric = self._calculate_metric_for_rna( topProteins)
            
            rnaShortestPath[ rna] = meanRNAShortestPath

            lionelMetrics[ rna] = lionelMetric

            ## calculate assortativity of graph, using current RNA scores

            # backup graph
            localGraph = self.graph.copy()

            # Add attribute to vertices in the graph, based on interaction with current RNA
            # note: if adding an attribute to one vertex, all other vertex will also have that attribute initialised as "None",
            # however, we do not have interactions for all proteins in PPI,
            # therefore the solution is to add a boolean, whether or not interaction occurs, based on a cutoff
            aboveCutoff = 0
            for scoreBin in self.rnaTargets[ rna]:
                for prot in self.rnaTargets[ rna][scoreBin]:
                    try:
                        vertexID = self.dictNames[ prot]
                        # if interaction exists, attribute value 1
                        if scoreBin > NetworkScoreAnalysis.SCORE_CUTOFF:
                            localGraph.vs[vertexID]["score"] = 1 #float( scoreBin)
                            aboveCutoff += 1
                        else:
                            localGraph.vs[vertexID]["score"] = 0    
                    except KeyError:
                        # not all proteins are in the PPI, this will be treated afterwards
                        continue
            # for the proteins in graph with no interaction, attribute value 0
            for v in localGraph.vs:
                if v["score"] == None:
                    v["score"] = 0

            assortativityMetric[ rna] = localGraph.assortativity("score", directed=False)
                    

            ## calculate metrics for each randomization

            if self.numberRandomizations > 0:

                if rna not in rnaShortestPathRandom:
                    rnaShortestPathRandom[ rna] = []
                if rna not in lionelMetricsRandom:
                    lionelMetricsRandom[ rna] = []
                if rna not in assortativityMetricRandom:
                    assortativityMetricRandom[ rna] = []
    
                for i in xrange( self.numberRandomizations):
                    
                    # retrieve a new set of top proteins, with the same degree
                    newProteinSet = self._get_sample_protein_degree( topProteins)
                    
                    meanRNAShortestPath, lionelMetric = self._calculate_metric_for_rna( newProteinSet)
                    
                    rnaShortestPathRandom[ rna].append( meanRNAShortestPath)
    
                    lionelMetricsRandom[ rna].append( lionelMetric)

                    ## randomise for assortativity metric
                    localGraph = self.graph.copy()

                    counter = 0                    
                    # randomize looping over graph
                    for v in random.sample( localGraph.vs, len( localGraph.vs)):
                        # add same number of interacting proteins as for the real case
                         
                        if counter < aboveCutoff:
                            v["score"] = 1
                            counter+=1
                        else:
                            v["score"] = 0


                    # testing
                    # value when all scores except one are 1: -8.10445996733e-06 (depends on which prot has the 0) # same is also true if I swap 1's and 0's
                    # value with aboveCutoff = 1000: -0.0690268123661
                    # value with aboveCutoff = 5000: -0.0923544580935
                    # value with aboveCutoff = 10000: -0.00756139885646
                    # 100 randomizations tend to have mean values of 0, regardless of aboveCutoff

                    assortativityMetricRandom[ rna].append( localGraph.assortativity("score", directed=False) )
    
                # calculate pvalue for lionel metric (based on randomization), they higher the better (count above)
                lionelMetricsPval[ rna] = self._empirical_pvalue( lionelMetricsRandom[ rna], lionelMetrics[ rna], 1)[0]
                # calculate pvalue for shortest path (based on randomization), they lower the better (count below)
                rnaShortestPathPval[ rna] = self._empirical_pvalue( rnaShortestPathRandom[ rna], rnaShortestPath[ rna], 0)[0]
                # calculate pvalue for assortativity (based on randomization), they higher the better (more correlated)
                assortativityMetricPval[ rna] = self._empirical_pvalue( assortativityMetricRandom[ rna], assortativityMetric[ rna], 1)[0]
     
                meanRandomLionelMetric = np.mean( lionelMetricsRandom[ rna])
                meanRandomRnaShortestPath = np.mean( rnaShortestPathRandom[ rna])
                meanRandomAssortativity = np.mean( assortativityMetricRandom[ rna])
                
                outFile.write( "%s\t%.2f\t%.2f\t%.1e\t%.2f\t%.2f\t%.1e\t%.3f\t%.3f\t%.1e\n" % ( rna, lionelMetrics[ rna], meanRandomLionelMetric, lionelMetricsPval[ rna], rnaShortestPath[ rna], meanRandomRnaShortestPath, rnaShortestPathPval[ rna], assortativityMetric[ rna], meanRandomAssortativity, assortativityMetricPval[ rna] ) )

            else:
                outFile.write( "%s\t%.2f\tNA\tNA\t%.2f\tNA\tNA\t%.3f\tNA\tNA\n" % ( rna, lionelMetrics[ rna], rnaShortestPath[ rna], assortativityMetric[ rna] ) )
                            
        outFile.close()
     
        #print (graph.average_path_length()) # take stime to run

        return lionelMetrics, lionelMetricsRandom, lionelMetricsPval, rnaShortestPath, rnaShortestPathRandom, rnaShortestPathPval


    # #
    # Internal function to actual perform the metric calculation for a given RNA
    # @param rna : ensembl ID of wanted RNA
    def _calculate_metric_for_rna(self, top_proteins):
        
        # get indexes of top proteins for wanted RNA        
        allIdx = [ self.dictNames[ prot] for prot in top_proteins]

        assert( len( allIdx) == self.topPartners)

        # stores mean shortest paths for this RNA, a value for each top protein
        meanShortestPaths = []

        # stores sum of lionel metric for this RNA
        lionelMetric = 0.0

        # for each top protein
        for idx in allIdx:
            
            # get list of indexes other than current
            otherIdx = [i for i in allIdx if i != idx ]
                            
            # calculate shortest path between current node against all others

            shortestPaths = []

#             # Approach 1: Calculating shortest path one vs one, saving result in a dictionary for reusal
#             for otherProt in otherIdx:
#                 tag = str( idx) + "|" + str( otherProt)
#                  
#                 if tag in self.calculatedShortestPaths:
#                     shortestPath = self.calculatedShortestPaths[ tag]
#                 else:
#                     shortestPath = self.graph.shortest_paths( idx, otherProt, mode = "OUT")[0]                
#                     # save result in dictionary. Shortest path is the same both ways in our type of network
#                     self.calculatedShortestPaths[ tag] = shortestPath
#                     self.calculatedShortestPaths[ str( otherProt) + "|" + str( idx)] = shortestPath
#                  
#                 shortestPaths.extend( shortestPath)
            
            # Approach 2: Calculating one protein vs all (much faster!)
            shortestPaths = self.graph.shortest_paths( idx, otherIdx, mode = "OUT")[0]
            
            #=======================================================================
            # Lionel metric
            #=======================================================================

            ## Shortest paths
            # distance 0 -> self interaction (should not happen in 'no self' network, also the shortest path is only calculated against other proteins)
            # distance 1 -> direct interaction
            # distance 2 -> one protein in the middle (first neighbours)
            # distance 3 -> two proteins in the middle (second neighbours)

            for pathLength in shortestPaths:
                
                if pathLength == 0:
                    raise RainetException( "NetworkScoreAnalysis.calculate_metrics : shortest path equals 0, error. %s" % ( idx ) )

                if pathLength in NetworkScoreAnalysis.DISTANCE_SCORE: 
                    #print pathLength, self.proteinGraphIDDict[ idx], self.proteinGraphIDDict[ otherIdx[ c]]
                    lionelMetric += NetworkScoreAnalysis.DISTANCE_SCORE[ pathLength]
                else:
                    lionelMetric += NetworkScoreAnalysis.DISTANCE_ABOVE_LIMIT_SCORE
           
            #=======================================================================
            # mean shortest path
            #=======================================================================
            
            # calculate mean shortest path for a node
            meanShortestPath = np.mean( shortestPaths)
            meanShortestPaths.append( meanShortestPath)

                    
        # calculate mean shortest path for current RNA (mean of means)
        meanRNAShortestPath = np.mean( meanShortestPaths)
                
        return meanRNAShortestPath, lionelMetric

    # #
    # Internal function to pick a random set of proteins with the same degree as input protein list.
    def _get_sample_protein_degree(self, list_of_proteins):
        
        degreeDict = self.degreeDict
        
        # make copy of proteins per degree in each randomization, so that I can remove already chosen proteins
        proteinsPerDegreeDictCopy = copy.deepcopy( self.proteinsPerDegreeDict)
        
        newProteinSet = set()

        # loop over provided list of proteins and remove them as a choice to be picked randomly
        for prot in list_of_proteins:
            degree = degreeDict[ prot]

            proteinsPerDegreeDictCopy[ degree].remove( prot)
            # if there is no more proteins with this degree, remove the degree as key in dictionary
            if len( proteinsPerDegreeDictCopy[ degree]) == 0:
                del proteinsPerDegreeDictCopy[ degree]

        # loop over provided list of proteins, pick and remove new set of proteins with same degree randomly
        for prot in list_of_proteins:
                        
            degree = degreeDict[ prot]

            # check if can pick protein with same degree, if not, sample protein from the closest possible degree
            if degree not in proteinsPerDegreeDictCopy:
                # change current degree to the closest available degree
                degree = self._get_new_degree( degree, proteinsPerDegreeDictCopy)
            elif len( proteinsPerDegreeDictCopy[ degree]) == 0:
                # change current degree to the closest available degree
                degree = self._get_new_degree( degree, proteinsPerDegreeDictCopy)

            # (meaning randomization is limited)
            if len( proteinsPerDegreeDictCopy[ degree]) < self.numberRandomizations:
                Logger.get_instance().warning( "NetworkScoreAnalysis._get_sample_protein_degree : less proteins with degree %s than number of randomizations." % ( degree ) )

            # pick one random protein from list of proteins with that degree
            randomProt = random.sample( proteinsPerDegreeDictCopy[ degree], 1)[0]

            # remove protein so that it is not picked up again in this randomization
            proteinsPerDegreeDictCopy[ degree].remove( randomProt)
            
            # if there is no more proteins with this degree, remove the degree as key in dictionary
            if len( proteinsPerDegreeDictCopy[ degree]) == 0:
                del proteinsPerDegreeDictCopy[ degree]
            
            newProteinSet.add( randomProt)
                   
        assert( len( list_of_proteins) == len( newProteinSet))
    
        # convert set into list
        newProteinSet = [prot for prot in newProteinSet]
    
        return newProteinSet


    # Function to get the closest possible available degree
    def _get_new_degree(self, degree, proteinsPerDegreeDictCopy):

        sortedDegrees = sorted( proteinsPerDegreeDictCopy.keys()) 
        
        # get list of all degrees except current degree
        otherDegrees = filter(lambda a: a != degree, sortedDegrees)

        # calculate distances of all degrees against ours and return the closest degree
        newDegree = min( otherDegrees, key=lambda x:abs( x-degree))

        Logger.get_instance().warning(  "NetworkScoreAnalysis._get_new_degree : No more proteins with degree %s. Sampling closest degree: %s." % ( degree, newDegree ) )

        return newDegree


    # #
    # Calculate proportion of random tests that are below observed value
    # Conservative approach: only counts observed below if < random value (not <=)
    def _empirical_pvalue(self, list_random, observed, aboveTail = 1):
        
        listRandomSignificant = np.sort( list_random)            

        below = 0 # count how many times observed value is below random values
        above = 0 # count how many times observed value is above random values
        for val in listRandomSignificant:
            if observed < val: # conservative
                below +=1
            if observed > val:
                above +=1

        if aboveTail:
            pval = 1.0 - ( above / float( len( listRandomSignificant)) )
            return pval, above
        else:
            pval = 1.0 - (below / float( len( listRandomSignificant)) )  
            return pval, below

    
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

        Timer.get_instance().step( "Calculate metrics.." )        

        self.calculate_metrics()
         
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
#         parser.add_argument('rainetDBFile', metavar='rainetDBFile', type=str,
#                              help='File with a RAINET Database to use for protein ID mapping.')
        parser.add_argument('topPartners', metavar='topPartners', type=int,
                             help='Number of top protein partners to consider for analysis, for each transcript.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Where to write output files.')
        # optional args
        parser.add_argument('--numberRandomizations', metavar='numberRandomizations', type=int, default = 1000,
                             help='Number of randomizations to be performed to calculate empirical p-value for each metric, for each transcript.')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        networkScoreAnalysis = NetworkScoreAnalysis( args.networkFile, args.catrapidFile, args.topPartners, args.outputFolder, args.numberRandomizations)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        networkScoreAnalysis.run( )

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

