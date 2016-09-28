import sys
import os
import argparse
import igraph

# import numpy as np
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
# 1)
# 2)
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

#         # make output folder
#         if not os.path.exists( self.outputFolder):
#             os.mkdir( self.outputFolder)


    # #
    # Use RAINET DB to retrieve Protein cross references
    def protein_cross_references(self):

        proteinIDMapping = "proteinIDMapping"

        DataManager.get_instance().perform_query( proteinIDMapping, "query( Protein.uniprotAC, Protein.uniprotID).all()") 

        # Convert query into a dictionary
        DataManager.get_instance().query_to_dict( proteinIDMapping, 1, 0)
        proteinIDMappingDict = DataManager.get_instance().get_data( proteinIDMapping) # key -> uniprotID (e.g. .._HUMAN), val -> uniprotAC (e.g. P35670)
        # formatting
        proteinIDMappingDict = { prot : proteinIDMappingDict[prot][0] for prot in proteinIDMappingDict}

        self.proteinIDMappingDict = proteinIDMappingDict

        return proteinIDMappingDict
   

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
                
        return graph, listOfNames, listOfTuples, dictNames
    
    
    # #
    # Function to calculate degree level for each protein in PPI network
    def calculate_protein_degree(self):

        degreeDict = {} # key -> uniprotAC, val -> degree of protein
        
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
                    
            try:
                proteinUniprotAC = self.proteinIDMappingDict[ proteinName]
            except KeyError:
                Logger.get_instance().warning( "NetworkScoreAnalysis.calculate_protein_degree : Could not find protein uniprotac for protein name %s. We will be converting here to uniprotAC" % ( proteinName ) )
                
                if "_HUMAN" in proteinName:
                    proteinUniprotAC = proteinName.split("_")[0]
                else:
                    proteinNotFound.add( proteinName)
                    continue
                        
            if proteinUniprotAC not in degreeDict:
                degreeDict[ proteinUniprotAC] = v.degree()
            else:
                raise RainetException( "NetworkScoreAnalysis.calculate_protein_degree : duplicate uniprotAC %s" % ( proteinUniprotAC ) )
 
        Logger.get_instance().info( "NetworkScoreAnalysis.calculate_protein_degree : Could not find proteins uniprotac for protein names %s. These were discarded." % ( proteinNotFound ) )
 
        self.degreeDict = degreeDict

        return degreeDict

    
    # #
    # Central function to run functions in order
    def run( self):

        Logger.get_instance().info( "NetworkScoreAnalysis.run: Starting..." )

        #===================================================================
        # Initialising datasets
        #===================================================================

        Timer.get_instance().step( "Reading network file.." )        
          
        self.read_network_file()


         
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
                             help='File with all vs all catrapid results. All interactions there will be processed (you may want to filter file in advance).')
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

