
import os
import shutil
import numpy
from scipy import stats

from sqlalchemy import or_, and_, distinct
from sqlalchemy.inspection import inspect    

from fr.tagc.rainet.core.execution.ExecutionStrategy import ExecutionStrategy
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.data.DataManager import DataManager
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil

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

from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util import Constants


# #
# This class define the Strategy to produce analysis on the given database and parameters
class EnrichmentAnalysisStrategy(ExecutionStrategy):

    #===============================================================================
    #
    # Enrichment Analysis strategy Constants
    #
    #===============================================================================

    # correspondance between wanted table and the base table
    ANNOTATION_TABLES_DICT = {"NetworkModule" : "ProteinNetworkModule", "ReactomePathway" : "ProteinReactomeAnnotation", "KEGGPathway" : "ProteinKEGGAnnotation"}



    #===================================================================
    # Data Manager object Keywords
    #===================================================================

    # Protein / RNA with interaction data
    PRI_PROT_KW = "interactingProteins" # Stores all Proteins in interactions
    PRI_RNA_KW = "interactingRNAs" # Stores RNAs in interactions

    # Proteins / RNAs with at least one interaction
    PRI_PROT_AT_LEAST_ONE_KW = "interactingProteinsAtLeastOne" # Stores all Proteins with at least one interaction
    PRI_RNA_AT_LEAST_ONE_KW = "interactingRNAsAtLeastOne" # Stores all RNAs with at least one interaction

    # Interactions
    PRI_KW = "interactions" # Stores interactions

    #===================================================================
    # Report files constants       
    #===================================================================

    PARAMETERS_LOG = "parameters.log"

    # Annotation report
    REPORT_PROT_PER_ANNOTATION = "prot_per_annotation.tsv"
    REPORT_ANNOTATION_PER_PROT = "annotation_per_prot.tsv"


    def __init__(self):  
        
        # Switch for writing of external report file      
        self.writeReportFile = 1


    # #
    # The Strategy execution method
    def execute(self):

        #===================================================================
        # Getting input arguments        
        #===================================================================
        
        self.DBPath = OptionManager.get_instance().get_option(OptionConstants.OPTION_DB_NAME)
        self.species = OptionManager.get_instance().get_option(OptionConstants.OPTION_SPECIES)
        self.outputFolder = OptionManager.get_instance().get_option(OptionConstants.OPTION_OUTPUT_FOLDER)
        self.annotationTable = OptionManager.get_instance().get_option(OptionConstants.OPTION_ANNOTATION_TABLE)

        # Variable that stores all arguments to appear in parameters log file
        self.arguments = {OptionConstants.OPTION_DB_NAME : self.DBPath,
                          OptionConstants.OPTION_SPECIES : self.species,
                          OptionConstants.OPTION_OUTPUT_FOLDER : self.outputFolder,
                          OptionConstants.OPTION_ANNOTATION_TABLE : self.annotationTable
                        }

        #===================================================================
        # Check input argument validity
        #===================================================================
        
        # Check if output folder path exists, create it if not
        if self.outputFolder != "" and len(self.outputFolder) > 0:
            FileUtils.initialise_output_folders(self.outputFolder)
        else:
            raise RainetException( "EnrichmentAnalysisStrategy.execute: Provided output folder is empty.")

        # check if annotation table name is consistent
        if self.annotationTable not in Constants.ANNOTATION_TABLES:
            raise RainetException( "EnrichmentAnalysisStrategy.execute: Provided annotation table name is not correct: " + self.annotationTable)
            

        #===================================================================
        # Initialisation
        #===================================================================

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.DBPath)
        self.sql_session = SQLManager.get_instance().get_session()
        self.db_engine = SQLManager.get_instance().get_engine()
                        
        self.analysis()
        
    # #
    # Central function to run analysis-related functions in order
    def analysis(self):

        Logger.get_instance().info( "EnrichmentAnalysisStrategy.analysis: Starting..." )

        Timer.get_instance().start_chrono()

        self.write_parameter_log()
          
        #===================================================================
        # Initialising datasets
        #===================================================================

        Timer.get_instance().step( "Initialising interaction datasets.." )        
          
        self.get_interaction_data()
          
        #===================================================================
        # Produce reports
        #===================================================================
    
        Timer.get_instance().step( "Producing annotation report.." )        
     
        self.annotation_report()
#   
#         #===================================================================
#         # Perform analysis
#         #===================================================================
#           
#         self.enrichement_analysis()
#           
#         if self.writeReportFile:
#             Timer.get_instance().step( "Writing report.." )
#             self.write_report()
# 
#         Timer.get_instance().stop_chrono( "Analysis Finished!")
 
    
    # #
    # Write output file with the parameters used
    def write_parameter_log(self):
                
        #===================================================================    
        # Write log of parameters used
        #=================================================================== 

        outHandler = FileUtils.open_text_w( self.outputFolder + "/" + EnrichmentAnalysisStrategy.PARAMETERS_LOG )
        
        Logger.get_instance().info( "\nARGUMENTS USED:" )

        outHandler.write( "Argument\tValue\n")
        for argName in sorted( self.arguments):
            argValue = self.arguments[ argName]
            Logger.get_instance().info( "%s:\t%s" % ( argName, argValue) ) 
            outHandler.write( "%s:\t%s\n" % ( argName, argValue) )
        outHandler.close()

    # #
    # Load items related to Protein-RNA interactions
    def get_interaction_data(self):

        # Get interacting proteins
        queryText = "query( InteractingProtein.uniprotAC ).all()"
        interactingProteins = eval('self.sql_session.' +  queryText)
        interactingProteins = { str( item[0]) for item in interactingProteins}

        Logger.get_instance().info( "get_interaction_data : Loaded %s interacting proteins. " % str( len( interactingProteins) ) )
        
        # Get interacting RNAs
        queryText = "query( InteractingRNA.transcriptID ).all()"
        interactingRNAs = eval('self.sql_session.' +  queryText)
        interactingRNAs = { str( item[0]) for item in interactingRNAs}
        
        Logger.get_instance().info( "get_interaction_data : Loaded %s interacting RNAs. " % str( len( interactingRNAs) ) )
        
        # Get interactions
        queryText = "query( ProteinRNAInteractionCatRAPID.transcriptID, ProteinRNAInteractionCatRAPID.proteinID, ProteinRNAInteractionCatRAPID.interactionScore ).all()"
        interactions = eval('self.sql_session.' +  queryText)

        Logger.get_instance().info( "get_interaction_data : Loaded %s interactions. " % str( len( interactions) ) )
        
        
#         # Get list of proteins with at least one interaction
#         setInterProts = set()
#         for interaction in interactions:
#             setInterProts.add( interaction.proteinID)
# 
#         Logger.get_instance().info( "get_interaction_data : proteins with at least one interaction. " % str( len( setInterProts) ) )
# 
#         # Get list of RNAs with at least one interaction
#         setInterRNAs = set()
#         for interaction in interactions:
#             setInterRNAs.add( interaction.transcriptID)
# 
#         Logger.get_instance().info( "get_interaction_data : proteins with at least one interaction. " % str( len( setInterRNAs) ) )       
        
        # store into accessible variables
        DataManager.get_instance().store_data(EnrichmentAnalysisStrategy.PRI_PROT_KW, interactingProteins)
        DataManager.get_instance().store_data(EnrichmentAnalysisStrategy.PRI_RNA_KW, interactingRNAs)
        DataManager.get_instance().store_data(EnrichmentAnalysisStrategy.PRI_KW, interactions)
        
#         DataManager.get_instance().store_data(EnrichmentAnalysisStrategy.PRI_PROT_KW, interactingProteins)
#         DataManager.get_instance().store_data(EnrichmentAnalysisStrategy.PRI_RNA_KW, interactingRNAs)


    # #
    # Retrieve statistics for the protein annotation data used.
    # Produce output files that will be used for a pdf report
    def annotation_report(self):
          
        #=================================================================== 
        # Retrieve interacting objects
        #=================================================================== 
 
        # proteins with interaction data
        allProteinsWithInteractionData = DataManager.get_instance().get_data( EnrichmentAnalysisStrategy.PRI_PROT_KW)
 
        # interactions after interaction filter
        interactions = DataManager.get_instance().get_data( EnrichmentAnalysisStrategy.PRI_KW)
  
        #===================================================================   
        # Retrieve annotations
        #=================================================================== 
        # Note: assuming there is two primary keys in annotation table, first pathway ID, then protein ID
        # Note: assuming there is only one primary key in base annotation table 
             
        # Define wanted tables
         
        # TODO: put annotation tables as system argument
         
        tableNameBase = self.annotationTable
        tableNameAnnotation = EnrichmentAnalysisStrategy.ANNOTATION_TABLES_DICT[ tableNameBase]
         
        # Get table primary key names 
        primaryKeys = eval( "inspect( " + tableNameAnnotation + ").primary_key") 
        # e.g. change from ProteinKEGGAnnotation.keggPathway_id to keggPathway_id
        primaryKeys = [ str(pk).replace( tableNameAnnotation, "")[ 1:] for pk in primaryKeys]
         
        # query table containing all pathway/module descriptions
        baseAnnotations = eval( "self.sql_session.query( " + tableNameBase +" ).all()")
        # query table containing all annotation mappings of pathway-protein
        proteinAnnotations = eval( "self.sql_session.query( " + tableNameAnnotation +" ).all()")
 
        #===================================================================   
        # Process annotations
        #=================================================================== 
         
        pathwayAnnotDict = {} # key -> pathway id, value -> list of proteins IDs
        protAnnotDict = {} # key -> protein id, value -> list of pathway IDs
 
        pathwayAnnotWithInteractionDataDict = {} # key -> pathway id, value -> list of proteins IDs with interaction data
        protAnnotWithInteractionDataDict = {} # key -> protein id (if with interaction data), value -> list of pathway IDs
 
        for annot in proteinAnnotations:
            pathID = str( eval( "annot." + primaryKeys[ 0] ) )
            protID = str( eval( "annot." + primaryKeys[ 1] ) )
 
            # populate dictionaries containing proteins with interaction data
            if protID in allProteinsWithInteractionData:
                if pathID not in pathwayAnnotWithInteractionDataDict:
                    pathwayAnnotWithInteractionDataDict[ pathID] = []
                pathwayAnnotWithInteractionDataDict[ pathID].append( protID)
                 
                if protID not in protAnnotWithInteractionDataDict:
                    protAnnotWithInteractionDataDict[ protID] = []
                if pathID not in protAnnotWithInteractionDataDict[ protID]:
                    protAnnotWithInteractionDataDict[ protID].append( pathID)
                else:
                    raise RainetException( "EnrichmentAnalysisStrategy.annotation_report: duplicate protein-annotation pair.")
 
            # populate dictionaries containing all proteins         
            if pathID not in pathwayAnnotDict:
                pathwayAnnotDict[ pathID] = []
            pathwayAnnotDict[ pathID].append( protID)
             
            if protID not in protAnnotDict:
                protAnnotDict[ protID] = []
            if pathID not in protAnnotDict[ protID]:
                protAnnotDict[ protID].append( pathID)
            else:
                raise RainetException( "EnrichmentAnalysisStrategy.annotation_report: duplicate protein-annotation pair.")
 
 
        # Note: the protein-pathway annotation table will only contain proteins and pathways which have any annotation
        # add here the pathways which did not have any protein annotation   
        # get primary key of base pathway table 
        baseAnnotationsPK = eval( "inspect( " + tableNameBase + ").primary_key") 
        # e.g. from ReactomePathway.reactomeID to reactomeID
        baseAnnotationsPK = str( baseAnnotationsPK[ 0]).replace( tableNameBase, "")[ 1:]
 
        for pathway in baseAnnotations:
            pathID = str( eval ( "pathway." + baseAnnotationsPK ) )
            if pathID not in pathwayAnnotDict:
                pathwayAnnotDict[ pathID] = []
            if pathID not in pathwayAnnotWithInteractionDataDict:
                pathwayAnnotWithInteractionDataDict[ pathID] = []
 
        assert ( len( pathwayAnnotDict) == len( pathwayAnnotWithInteractionDataDict) )
 
#         print ( len( pathwayAnnotDict)) # 300. Put this in unittest?
#         print ( len( protAnnotDict))          
#         print pathwayAnnotDict
#         print pathwayAnnotWithInteractionDataDict
#         print protAnnotDict
#         print protAnnotWithInteractionDataDict
#          
        #===================================================================   
        # File with proteins per pathway
        #=================================================================== 
        # Note: not all proteins with interaction have annotation, and vice versa
 
        outHandler = FileUtils.open_text_w( self.outputFolder + "/" + EnrichmentAnalysisStrategy.REPORT_PROT_PER_ANNOTATION )
  
        # Write header
        outHandler.write("annotationID\ttotal_prot_in_pathway\ttotal_prot_with_interaction_data\n") 
  
        for pathID in pathwayAnnotDict:
 
            annot = pathwayAnnotDict[ pathID]
            annotWithInteractionData = pathwayAnnotWithInteractionDataDict[ pathID]
             
            outHandler.write( "%s\t%s\t%s\n" % ( pathID, len( annot), len( annotWithInteractionData)) )
          
        outHandler.close()
 
#         #===================================================================   
#         # File with pathway per protein
#         #===================================================================          
#  
#         outHandler = FileUtils.open_text_w( self.outputFolder + "/" + EnrichmentAnalysisStrategy.REPORT_ANNOTATION_PER_PROT )
#   
#         # Write header
#         outHandler.write("uniprotAC\twith_interaction_data\ttotal_annotations\tlist_of_annotations\n") 
#   
#         for protID in protAnnotDict:
#              
#             withInteraction = 0
#             if protID in interProtObjects:
#                 withInteraction = 1
#  
#             annot = protAnnotDict[ protID]
#                          
#             outHandler.write( "%s\t%s\t%s\t%s\n" % ( protID, withInteraction, len( annot), ",".join( annot) ) )
#           
#         outHandler.close()
 
 
#         #### TESTING / EXPLORATORY ####
#   
#         # writing file with one line per RNA - module pair, if interactions are existing
#   
#         outFile = open("/home/diogo/testing/networkModules/networkModules.tsv","w")
#   
#         outFile.write("transcriptID\tannotID\tnumber_interactions\tnumber_possible_interactions\ttotal_interacting_proteins\tlist_proteinIDs\n")
#   
#         # For each RNA, store all proteins it interacts with and their annotations
#         rnaInteractions = {} # Key -> transcript ID, value -> dict; key -> pathway ID, value -> list of prot IDs (after filtering)
#         totalRNAInteractionsDict = {} # Key -> transcript ID, value -> list of prot IDs (after filtering)       
#         for inter in interactions:
#             txID = str( inter.transcriptID)
#             protID = str( inter.proteinID)
#   
#             if txID not in rnaInteractions:
#                 rnaInteractions[ txID] = {}
#                 totalRNAInteractionsDict[ txID] = []
#             totalRNAInteractionsDict[ txID].append( protID)
#               
#             # only store info of proteins that have annotation information
#             if protID in protAnnotDict: 
#                 for annot in protAnnotDict[ protID]:
#                     if annot not in rnaInteractions[ txID]:
#                         rnaInteractions[ txID][ annot] = []
#                     rnaInteractions[ txID][ annot].append( protID)
#   
#   
#         print "RNAs with interactions:", len( rnaInteractions)
#   
#         for rnaID in sorted( rnaInteractions):
# #            totAnnotatedInteractions = sum( [ len( rnaInteractions[ rnaID][ annotID]) for annotID in rnaInteractions[ rnaID] ] ) 
#             totRNAInteractions = len( totalRNAInteractionsDict[ rnaID]) 
#              
#             for annotID in sorted( rnaInteractions[ rnaID]):
#                  
# #                 # TODO: add minimum number of RBPs cutoff as argument or something???                
# #                 if len( pathwayAnnotDict[ annotID]) < 4: # TO REMOVE
# #                     continue
#                  
#                 protList = rnaInteractions[ rnaID][ annotID]
#                 possibleProtList = pathwayAnnotWithInteractionDataDict[ annotID]
#                   
#                 assert ( len( protList) <= len( possibleProtList) )
#   
#                 # note: test for each RNA - annotation pair. 
#                 x = len( protList) # white balls drawn ( proteins with current annotation in filtered interactions)
#                 m = len( possibleProtList) # total white balls ( proteins with the current annotation)
#                 n = len( interProtObjects) - m # total black balls (all the proteins used in catRAPID)
#                 k = totRNAInteractions # total number of draws ( proteins in filtered interactions)
#   
#                 assert ( m + n >= k) # number of draws should be less than total number of balls
#   
#                 hyperResult = self.hypergeometric_test(x, m, n, k)
#  
#                 #print ( rnaID, annotID, len( protList), len( possibleProtList), totRNAInteractions, ",".join( sorted(protList) ) ) 
#                   
#                 outFile.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( rnaID, annotID, len( protList), len( possibleProtList), totRNAInteractions, ",".join( sorted(protList) ), hyperResult ) )
#  
#         outFile.close()

        
    def enrichement_analysis(self):
        pass


    # #
    # Run hypergeometric test
    def hypergeometric_test(self, x, m, n, k):

        # x, q vector of quantiles representing the number of white balls drawn
        # without replacement from an urn which contains both black and white
        # balls.
        # 
        # m the number of white balls in the urn.
        # 
        # n the number of black balls in the urn.
        # 
        # k the number of balls drawn from the urn.

        # stats.hypergeome.sf gives the same result as R phyper
        testResult = stats.hypergeom.sf( x, m+n, m, k)

        # print ("x: %i\tm: %i\tn: %i\tk: %i" % (x,m,n,k) )
        # print ("Hypergeometric_test p-value:\t%.3f\n" % testResult )

        return testResult


    # #
    # Run Rscript to produce Sweave file and consequent pdf report, using the data written by this script
    def write_report(self):
        pass

