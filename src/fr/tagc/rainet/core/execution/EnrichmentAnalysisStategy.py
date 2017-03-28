
import os
import shutil
import numpy
import random
from scipy import stats
# from statsmodels.sandbox.stats.multicomp import multipletests

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
from fr.tagc.rainet.core.data.BioplexCluster import BioplexCluster
from fr.tagc.rainet.core.data.ProteinBioplexAnnotation import ProteinBioplexAnnotation
from fr.tagc.rainet.core.data.WanCluster import WanCluster
from fr.tagc.rainet.core.data.ProteinWanAnnotation import ProteinWanAnnotation
from fr.tagc.rainet.core.data.CorumCluster import CorumCluster
from fr.tagc.rainet.core.data.ProteinCorumAnnotation import ProteinCorumAnnotation
from fr.tagc.rainet.core.data.CustomCluster import CustomCluster
from fr.tagc.rainet.core.data.ProteinCustomAnnotation import ProteinCustomAnnotation


from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util import Constants
from statsmodels.stats.multitest import multipletests


# #
# This class define the Strategy to produce analysis on the given database and parameters
class EnrichmentAnalysisStrategy(ExecutionStrategy):

    #===============================================================================
    # Enrichment Analysis strategy Constants
    #===============================================================================

    # correspondance between the base table and associated data
    ANNOTATION_TABLES_DICT = {"NetworkModule" : "ProteinNetworkModule",
                              "ReactomePathway" : "ProteinReactomeAnnotation",
                              "KEGGPathway" : "ProteinKEGGAnnotation",
                              "BioplexCluster" : "ProteinBioplexAnnotation",
                              "WanCluster" : "ProteinWanAnnotation",
                              "CorumCluster" : "ProteinCorumAnnotation",
                              "CustomCluster" : "ProteinCustomAnnotation"}

    # significance value 
    SIGN_VALUE_TEST = 0.05  # the alpha value to call enrichment for each individual hypergeomtric test
    SIGN_VALUE_AGAINST_RANDOM_CONTROL = 0.01  # the alpha value for calling a observed different than random control
        
    SIGN_COLUMN = 8  # column from testsCorrected with result significance tag. 0-based
    WARNING_COLUMN = 5  # column from testsCorrected with warning/skipTest tag. 0-based   

    #===================================================================
    # Data Manager object Keywords
    #===================================================================

    # Protein / RNA with interaction data
    PRI_PROT_KW = "interactingProteins"  # Stores all Proteins in interactions
    PRI_RNA_KW = "interactingRNAs"  # Stores RNAs in interactions

    # Proteins / RNAs with at least one interaction
    PRI_PROT_AT_LEAST_ONE_KW = "interactingProteinsAtLeastOne"  # Stores all Proteins with at least one interaction
    PRI_RNA_AT_LEAST_ONE_KW = "interactingRNAsAtLeastOne"  # Stores all RNAs with at least one interaction

    # Interactions
    PRI_KW = "interactions"  # Stores interactions

    #===================================================================
    # Report files constants       
    #===================================================================

    PARAMETERS_LOG = "parameters.log"

    # Annotation report
    REPORT_PROT_PER_ANNOTATION = "prot_per_annotation.tsv"
    REPORT_ANNOTATION_PER_PROT = "annotation_per_prot.tsv"
    REPORT_ENRICHMENT = "enrichment_results.tsv"
    REPORT_ENRICHMENT_PER_RNA = "enrichment_per_rna.tsv"


    def __init__(self):  
        
        # Switch for writing of external report file      
        self.writeReportFile = 1

        # Container for results of hypergeom tests so that they don't have to be repeated. 
        self.testContainer = {}  # Key -> test parameters, val -> pval result

        # counter for total amount of tests
        self.countTotalTests = 0

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
        self.minimumProteinAnnotation = OptionManager.get_instance().get_option(OptionConstants.OPTION_MINIMUM_PROTEIN_ANNOTATION)
        self.minimumProteinInteraction = OptionManager.get_instance().get_option(OptionConstants.OPTION_MINIMUM_PROTEIN_INTERACTION)
        self.numberRandomizations = OptionManager.get_instance().get_option(OptionConstants.OPTION_NUMBER_RANDOMIZATIONS)
        self.expressionWarning = OptionManager.get_instance().get_option(OptionConstants.OPTION_EXPRESSION_WARNING)
        self.minimumExpression = OptionManager.get_instance().get_option(OptionConstants.OPTION_MINIMUM_EXPRESSION)
        self.lowerTail = OptionManager.get_instance().get_option(OptionConstants.OPTION_LOWER_TAIL)

        # Variable that stores all arguments to appear in parameters log file
        self.arguments = {OptionConstants.OPTION_DB_NAME : self.DBPath,
                          OptionConstants.OPTION_SPECIES : self.species,
                          OptionConstants.OPTION_OUTPUT_FOLDER : self.outputFolder,
                          OptionConstants.OPTION_ANNOTATION_TABLE : self.annotationTable,
                          OptionConstants.OPTION_MINIMUM_PROTEIN_ANNOTATION : self.minimumProteinAnnotation,
                          OptionConstants.OPTION_MINIMUM_PROTEIN_INTERACTION : self.minimumProteinInteraction,
                          OptionConstants.OPTION_NUMBER_RANDOMIZATIONS : self.numberRandomizations,
                          OptionConstants.OPTION_EXPRESSION_WARNING : self.expressionWarning,
                          OptionConstants.OPTION_MINIMUM_EXPRESSION : self.minimumExpression,
                          OptionConstants.OPTION_LOWER_TAIL : self.lowerTail                          
                        }

        #===================================================================
        # Check input argument validity
        #===================================================================
        
        # Check if output folder path exists, create it if not
        if self.outputFolder != "" and len(self.outputFolder) > 0:
            FileUtils.initialise_output_folders(self.outputFolder)
        else:
            raise RainetException("EnrichmentAnalysisStrategy.execute: Provided output folder is empty.")

        # check if annotation table name is consistent
        if self.annotationTable not in Constants.ANNOTATION_TABLES:
            raise RainetException("EnrichmentAnalysisStrategy.execute: Provided annotation table name is not correct: " + self.annotationTable + "Must be one of : " + str(Constants.ANNOTATION_TABLES))
            
        # Check if minimum protein annotation value is consistent
        try:
            self.minimumProteinAnnotation = int(self.minimumProteinAnnotation)
        except:
            raise RainetException("EnrichmentAnalysisStrategy.execute: Provided minimum protein annotation value is not correct (must be integer): " + str( self.minimumProteinAnnotation))
            
        # Check if minimum protein annotation value is consistent
        try:
            self.minimumProteinInteraction = int(self.minimumProteinInteraction)
        except:
            raise RainetException("EnrichmentAnalysisStrategy.execute: Provided minimum protein interactions value is not correct (must be integer): " + str( self.minimumProteinInteraction))

        # Check if minimum protein annotation value is consistent
        try:
            self.numberRandomizations = int(self.numberRandomizations)
        except:
            raise RainetException("EnrichmentAnalysisStrategy.execute: Provided number randomizations value is not correct (must be integer): " + str( self.numberRandomizations))

        # Check if expression warning is consistent
        if self.expressionWarning != OptionConstants.DEFAULT_EXPRESSION_WARNING:            
            try:
                self.expressionWarning = float(self.expressionWarning)
                if self.expressionWarning < 0.0 or self.expressionWarning > 1.0:
                    raise RainetException("EnrichmentAnalysisStrategy.execute: Provided expression warning value is not correct (must be float between 0.0 and 1.0): " + str(self.expressionWarning))
            except:
                raise RainetException("EnrichmentAnalysisStrategy.execute: Provided expression warning value is not correct (must be float between 0.0 and 1.0): " + str(self.expressionWarning))

        # Check if minimum protein annotation value is consistent
        if self.minimumExpression != OptionConstants.DEFAULT_MINIMUM_EXPRESSION and self.expressionWarning == OptionConstants.DEFAULT_EXPRESSION_WARNING:            
            raise RainetException("EnrichmentAnalysisStrategy.execute: Provided minimum RPKM expression value only works if --expressionWarning is active.")
        try:
            self.minimumExpression = float(self.minimumExpression)
        except:
            raise RainetException("EnrichmentAnalysisStrategy.execute: Provided minimum RPKM expression value is not correct (must be float): " + str( self.minimumExpression))

        # Check if minimum protein annotation value is consistent
        try:
            self.lowerTail = int(self.lowerTail)
        except:
            raise RainetException("EnrichmentAnalysisStrategy.execute: Provided lower tail option is not correct (must be integer): " + str( self.lowerTail))

        #===================================================================
        # Initialisation
        #===================================================================

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.DBPath)
        self.sql_session = SQLManager.get_instance().get_session()
        self.db_engine = SQLManager.get_instance().get_engine()
         
        ##TODO: UNCOMMENT THIS!!               
        self.analysis()
        
    # #
    # Central function to run analysis-related functions in order
    def analysis(self):

        Logger.get_instance().info("EnrichmentAnalysisStrategy.analysis: Starting...")

        Timer.get_instance().start_chrono()

        self.write_parameter_log()
          
        #===================================================================
        # Initialising datasets
        #===================================================================

        Timer.get_instance().step("Initialising interaction datasets..")        
          
        self.get_interaction_data()
          
        #===================================================================
        # Produce reports
        #===================================================================
    
        Timer.get_instance().step("Producing annotation report..")        
     
        self.annotation_report()
           

        #===================================================================
        # Load expression data (if needed)
        #===================================================================
           
        if self.expressionWarning != OptionConstants.DEFAULT_EXPRESSION_WARNING:
            
            Timer.get_instance().step("Producing annotation report..")        
            
            self.retrieve_expression()
   
        #===================================================================
        # Perform analysis
        #===================================================================

        Timer.get_instance().step("Producing enrichment report..")        
           
        self.enrichement_analysis()
#           
#         if self.writeReportFile:
#             Timer.get_instance().step( "Writing report.." )
#             self.write_report()
# 
        # Container for results of hypergeom tests so that they don't have to be repeated. 
        Logger.get_instance().info("EnrichmentAnalysisStrategy.analysis : Number of unique tests performed: %s" % (len(self.testContainer)))
        Logger.get_instance().info("EnrichmentAnalysisStrategy.analysis : Number of tests performed: %s" % (self.countTotalTests))
        
        Timer.get_instance().stop_chrono("Analysis Finished!")
 
    
    # #
    # Write output file with the parameters used
    def write_parameter_log(self):
                
        #===================================================================    
        # Write log of parameters used
        #=================================================================== 

        outHandler = FileUtils.open_text_w(self.outputFolder + "/" + EnrichmentAnalysisStrategy.PARAMETERS_LOG)
        
        Logger.get_instance().info("\nARGUMENTS USED:")

        outHandler.write("Argument\tValue\n")
        for argName in sorted(self.arguments):
            argValue = self.arguments[ argName]
            Logger.get_instance().info("%s:\t%s" % (argName, argValue)) 
            outHandler.write("%s:\t%s\n" % (argName, argValue))
        outHandler.close()

    # #
    # Load items related to Protein-RNA interactions
    def get_interaction_data(self):

        # Get interacting proteins
        queryText = "query( InteractingProtein.uniprotAC ).all()"
        interactingProteins = eval('self.sql_session.' + queryText)
        interactingProteins = { str(item[0]) for item in interactingProteins}

        Logger.get_instance().info("get_interaction_data : Loaded %s interacting proteins. " % str(len(interactingProteins)))
        
        # Get interacting RNAs
        queryText = "query( InteractingRNA.transcriptID ).all()"
        interactingRNAs = eval('self.sql_session.' + queryText)
        interactingRNAs = { str(item[0]) for item in interactingRNAs}
        
        Logger.get_instance().info("get_interaction_data : Loaded %s interacting RNAs. " % str(len(interactingRNAs)))
        
        # Get interactions
        queryText = "query( ProteinRNAInteractionCatRAPID.transcriptID, ProteinRNAInteractionCatRAPID.proteinID, ProteinRNAInteractionCatRAPID.interactionScore ).all()"
        interactions = eval('self.sql_session.' + queryText)

        Logger.get_instance().info("get_interaction_data : Loaded %s interactions. " % str(len(interactions)))
        
        
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
        # Retrieve interacting data
        #=================================================================== 
 
        # proteins with interaction data, based on interactingProteins table of RAINET DB
        allProteinsWithInteractionData = DataManager.get_instance().get_data(EnrichmentAnalysisStrategy.PRI_PROT_KW)
  
        #===================================================================   
        # Retrieve annotations
        #=================================================================== 
        # Note: assuming there is two primary keys in annotation table, first pathway ID, then protein ID
        # Note: assuming there is only one primary key in base annotation table 
             
        # Define wanted tables         
        tableNameBase = self.annotationTable
        tableNameAnnotation = EnrichmentAnalysisStrategy.ANNOTATION_TABLES_DICT[ tableNameBase]
         
        # Get table primary key names 
        primaryKeys = eval("inspect( " + tableNameAnnotation + ").primary_key") 
        # e.g. change from ProteinKEGGAnnotation.keggPathway_id to keggPathway_id
        primaryKeys = [ str(pk).replace(tableNameAnnotation, "")[ 1:] for pk in primaryKeys]
         
        # query table containing all pathway/module descriptions
        baseAnnotations = eval("self.sql_session.query( " + tableNameBase + " ).all()")
        # query table containing all annotation mappings of pathway-protein
        proteinAnnotations = eval("self.sql_session.query( " + tableNameAnnotation + " ).all()")
  
        #===================================================================   
        # Process annotations
        #=================================================================== 
         
        pathwayAnnotDict = {}  # key -> pathway id, value -> list of proteins IDs
        protAnnotDict = {}  # key -> protein id, value -> list of pathway IDs
 
        pathwayAnnotWithInteractionDataDict = {}  # #FOR REPORTING PURPOSES key -> pathway id, value -> list of proteins IDs with interaction data
        protAnnotWithInteractionDataDict = {}  # key -> protein id (if with interaction data), value -> list of pathway IDs
        annotWithInteractionDict = {}  # #FOR ANALYSIS PURPOSES. key -> pathway id (if with interaction data), value -> list of proteins IDs with interaction data

        for annot in proteinAnnotations:
            pathID = str(eval("annot." + primaryKeys[ 0]))
            protID = str(eval("annot." + primaryKeys[ 1]))
 
            # populate dictionaries containing proteins with interaction data
            if protID in allProteinsWithInteractionData:
                
                # initialise if interaction only
                if pathID not in annotWithInteractionDict:
                    annotWithInteractionDict[ pathID] = []
                annotWithInteractionDict[ pathID].append(protID)

                if pathID not in pathwayAnnotWithInteractionDataDict:
                    pathwayAnnotWithInteractionDataDict[ pathID] = []
                pathwayAnnotWithInteractionDataDict[ pathID].append(protID)
                 
                if protID not in protAnnotWithInteractionDataDict:
                    protAnnotWithInteractionDataDict[ protID] = []

                if pathID not in protAnnotWithInteractionDataDict[ protID]:
                    protAnnotWithInteractionDataDict[ protID].append(pathID)
                else:
                    raise RainetException("EnrichmentAnalysisStrategy.annotation_report: duplicate protein-annotation pair.")
 
            # populate dictionaries containing all proteins
            # initialise regardless of interaction
            if pathID not in pathwayAnnotDict:
                pathwayAnnotDict[ pathID] = []
            pathwayAnnotDict[ pathID].append(protID)
             
            if protID not in protAnnotDict:
                protAnnotDict[ protID] = []
            if pathID not in protAnnotDict[ protID]:
                protAnnotDict[ protID].append(pathID)
            else:
                raise RainetException("EnrichmentAnalysisStrategy.annotation_report: duplicate protein-annotation pair.")
  
        # Note: the protein-pathway annotation table will only contain proteins and pathways which have any annotation
        # add here the pathways which did not have any protein annotation   
        # get primary key of base pathway table 
        baseAnnotationsPK = eval("inspect( " + tableNameBase + ").primary_key") 
        # e.g. from ReactomePathway.reactomeID to reactomeID
        baseAnnotationsPK = str(baseAnnotationsPK[ 0]).replace(tableNameBase, "")[ 1:]
  
        for pathway in baseAnnotations:
            pathID = str(eval ("pathway." + baseAnnotationsPK))
            if pathID not in pathwayAnnotDict:
                pathwayAnnotDict[ pathID] = []
            if pathID not in pathwayAnnotWithInteractionDataDict:
                pathwayAnnotWithInteractionDataDict[ pathID] = []
 
        assert (len(pathwayAnnotDict) == len(pathwayAnnotWithInteractionDataDict))

        #===================================================================   
        # File with proteins per pathway
        #=================================================================== 
        # Note: not all proteins with interaction have annotation, and vice versa
 
        outHandler = FileUtils.open_text_w(self.outputFolder + "/" + EnrichmentAnalysisStrategy.REPORT_PROT_PER_ANNOTATION)
  
        # Write header
        outHandler.write("annotationID\ttotal_prot_in_pathway\ttotal_prot_with_interaction_data\n") 
  
        summ = 0
        summWithInteractionData = 0
        for pathID in pathwayAnnotDict:
 
            annot = set(pathwayAnnotDict[ pathID])
            annotWithInteractionData = set(pathwayAnnotWithInteractionDataDict[ pathID]) 
            
            summ += len(annot)
            summWithInteractionData += len(annotWithInteractionData)
            
            assert(len(annotWithInteractionData) == len(annot.intersection(annotWithInteractionData))), "annotations with interaction need to be present in all-annotations dictionary"
             
            outHandler.write("%s\t%s\t%s\n" % (pathID, len(annot), len(annotWithInteractionData)))

        outHandler.write("# Annotation coverage: %.2f%%\n" % (summWithInteractionData * 100.0 / summ))
          
        outHandler.close()

        Logger.get_instance().info("EnrichmentAnalysisStrategy.annotation_report: Number annotations: %s " % str(len(pathwayAnnotDict)))
        Logger.get_instance().info("EnrichmentAnalysisStrategy.annotation_report: Number annotations with interactions: %s " % str(len(annotWithInteractionDict)))
 
        #===================================================================   
        # File with pathway per protein
        #===================================================================          
  
        outHandler = FileUtils.open_text_w(self.outputFolder + "/" + EnrichmentAnalysisStrategy.REPORT_ANNOTATION_PER_PROT)
   
        # Write header
        outHandler.write("uniprotAC\twith_interaction_data\ttotal_annotations\tlist_of_annotations\n") 
      
        sumWithInteractionData = 0
        for protID in protAnnotDict:
              
            withInteraction = 0
            if protID in allProteinsWithInteractionData:
                withInteraction = 1
                sumWithInteractionData += 1
  
            annot = protAnnotDict[ protID]
                          
            outHandler.write("%s\t%s\t%s\t%s\n" % (protID, withInteraction, len(annot), ",".join(annot)))

        Logger.get_instance().info("EnrichmentAnalysisStrategy.annotation_report: # Proteins with at least one annotation: %i. # with interaction data: %i" % (len(protAnnotDict), sumWithInteractionData))
        outHandler.write("# # Proteins with at least one annotation: %i. # with interaction data: %i\n" % (len(protAnnotDict), sumWithInteractionData))
           
        outHandler.close()
         
        #===================================================================   
        # store processed data
        #===================================================================   
        self.protAnnotDict = protAnnotDict
        self.annotWithInteractionDict = annotWithInteractionDict
        self.allProteinsWithInteractionDataLen = len(allProteinsWithInteractionData)
        
        # proteins with interaction and with annotation # for expression filtering
        self.proteinWithAnnotationWithInteraction = { prot for pathID in annotWithInteractionDict for prot in annotWithInteractionDict[ pathID]}

        assert sumWithInteractionData == len(self.proteinWithAnnotationWithInteraction)

        # Logger.get_instance().info( "EnrichmentAnalysisStrategy.annotation_report: Proteins with interactions, based on interactingProteins table: %s " % str( self.allProteinsWithInteractionDataLen ) )
        Logger.get_instance().info("EnrichmentAnalysisStrategy.annotation_report: Proteins with interactions and with annotations: %s " % str(len(self.proteinWithAnnotationWithInteraction)))

 
    # #
    # Perform enrichment analysis
    # Write main output files
    def enrichement_analysis(self):

        #=================================================================== 
        # Retrieve interactions
        #=================================================================== 
         
        # interactions after interaction filter
        interactions = DataManager.get_instance().get_data(EnrichmentAnalysisStrategy.PRI_KW)

        # For each RNA, store all proteins it interacts with and their annotations
        rnaInteractions = {}  # Key -> transcript ID, value -> dict; key -> pathway ID, value -> list of prot IDs (after filtering)
        
        setInteractingProteins = set()  # stores proteins with at least one interaction
        
        for inter in interactions:
            txID = str(inter.transcriptID)
            protID = str(inter.proteinID)
            
            setInteractingProteins.add(protID)
                               
            # only store info of proteins that have annotation information
            if protID in self.protAnnotDict: 

                # only initialise RNA in dictionary if there is at least one protein with annotation   
                if txID not in rnaInteractions:
                    rnaInteractions[ txID] = {}

                for annot in self.protAnnotDict[ protID]:
                    if annot not in rnaInteractions[ txID]:
                        rnaInteractions[ txID][ annot] = []
                    rnaInteractions[ txID][ annot].append(protID)
      
        Logger.get_instance().info("EnrichmentAnalysisStrategy.enrichement_analysis: RNAs with interactions: %s " % str(len(rnaInteractions)))
        Logger.get_instance().info("EnrichmentAnalysisStrategy.enrichement_analysis: Proteins with at least one interaction: %s " % str(len(setInteractingProteins)))
        Logger.get_instance().info("EnrichmentAnalysisStrategy.enrichement_analysis: Proteins with annotation: %s " % str(len(self.protAnnotDict)))

        self.rnaInteractionsLength = len(rnaInteractions)

        #===================================================================   
        # Set protein background
        #===================================================================   
        # Background (search space) changes depending on annotation dataset being used
        # it is the intersection of proteins with at least one interaction and proteins with annotation

        backgroundProteins = setInteractingProteins.intersection(set(self.protAnnotDict.keys()))
        self.backgroundProteins = backgroundProteins

        Logger.get_instance().info("EnrichmentAnalysisStrategy.enrichement_analysis: Protein background: %s " % str(len(backgroundProteins)))

        # Note: difference between backgroundProteins and 'self.proteinWithAnnotationWithInteraction' from annotation_report, 
        # is that in annotatio_report the set of interacting proteins is retrieved from interactingProteins table, 
        # whereas here it is retrieved straight from the interaction table, depending on the interaction filtering applied,
        # not all proteins in interactingProteins table will contain at least an interacting in the interaction table.

        # delete object to save memory
        self.protAnnotDictLen = len(self.protAnnotDict)
        del self.protAnnotDict

        #===================================================================    
        #
        # Make enrichment tests
        # 
        #===================================================================          
  
        #===================================================================   
        # Initialise file with enrichment test per RNA-module pair
        #===================================================================          
        # one line per RNA-module (enrichment test)
        outHandler = FileUtils.open_text_w(self.outputFolder + "/" + EnrichmentAnalysisStrategy.REPORT_ENRICHMENT)
        outHandler.write("transcriptID\tannotID\tnumber_observed_interactions\tnumber_possible_interactions\tpercent_interacting_proteins\twarning\tpval\tcorrected_pval\tsign_corrected\n")

        #===================================================================   
        # Initialise file with stats per RNA
        #===================================================================          
        # one line per RNA, observed significant tests vs random significant tests
        outHandlerStats = FileUtils.open_text_w(self.outputFolder + "/" + EnrichmentAnalysisStrategy.REPORT_ENRICHMENT_PER_RNA)
        outHandlerStats.write("transcriptID\tn_sign_tests_no_warning\tavg_n_sign_random_no_warning\tn_times_above_random\tempiricalPvalue\tsignificant\n")

        #===================================================================   
        # Annotation randomization
        #===================================================================          

        # Calculate number of possible permutations
        listOfProteins = [ prot for annot in self.annotWithInteractionDict for prot in self.annotWithInteractionDict[ annot]]
        avgPotSize = len(listOfProteins) / float(len(self.annotWithInteractionDict))
        nPermutations = avgPotSize * len(listOfProteins)
        Logger.get_instance().info("EnrichmentAnalysisStrategy.enrichement_analysis : Number of possible permutations: %i" % (nPermutations))
        if nPermutations < self.numberRandomizations:
            Logger.get_instance().warning("EnrichmentAnalysisStrategy.enrichement_analysis : Number of possible permutations is less than number of randomizations.")           
       
        # Get list of proteins in a orderly manner
        listOfProteins = [ prot for annot in sorted(self.annotWithInteractionDict) for prot in sorted(self.annotWithInteractionDict[ annot]) ]
       
        # shuffle annotation tags of proteins
        randomAnnotDicts = [ self.randomize_proteins(self.annotWithInteractionDict, listOfProteins) for i in xrange(0, self.numberRandomizations)]

        #===================================================================   
        #===================================================================   
        # Loop RNA to perform tests per RNA
        #===================================================================          
        #===================================================================   
        rnaCounter = 0
        # for each RNA with any interaction
        for rnaID in sorted(rnaInteractions):
                        
            rnaCounter += 1
            if rnaCounter % 100 == 0:
                Timer.get_instance().step("EnrichmentAnalysisStrategy.enrichement_analysis : processed %s RNAs.." % str(rnaCounter))        

            # retrieve total number of interactions with annotated proteins for this RNA
            totalRNAInteractions = { prot for annotID in rnaInteractions[ rnaID] for prot in rnaInteractions[ rnaID][ annotID]}

            #===================================================================   
            # Run real test            
            #=================================================================== 
            
            if self.expressionWarning != OptionConstants.DEFAULT_EXPRESSION_WARNING:            
                if rnaID not in self.expressionDict:
                    Logger.get_instance().warning("enrichement_analysis : skipping rna due to lack of expression data: %s" % (rnaID))
                    continue
              
            testsCorrected = self.run_rna_vs_annotations(rnaID, self.annotWithInteractionDict, totalRNAInteractions)

            assert len(testsCorrected) == len(self.annotWithInteractionDict)

            #===================================================================          
            #===================================================================   
            # File with enrichment test per RNA-module pair
            #===================================================================          
            #===================================================================          

            for test in testsCorrected:
                outHandler.write("\t".join(test) + "\n")

            #===================================================================   
            # Randomization tests
            #===================================================================          

            listRandomSignificants = numpy.empty(self.numberRandomizations, object)
            listRandomSignificantsNoWarning = numpy.empty(self.numberRandomizations, object)
            for i in xrange(0, self.numberRandomizations):               
                randomTestsCorrected = self.run_rna_vs_annotations(rnaID, randomAnnotDicts[ i], totalRNAInteractions)           
                listRandomSignificants[ i], listRandomSignificantsNoWarning[ i] = self.count_sign_tests(randomTestsCorrected)

            # using just Significant no warning
            avgSignRandomNoWarning = numpy.mean(listRandomSignificantsNoWarning)

            #===================================================================          
            #===================================================================   
            # File with stats per RNA
            #===================================================================          
            #===================================================================          
            # count number of significant tests AFTER p-value correction
            
            # stats for the real/observed data
            countSignificant, countSignificantNoWarning = self.count_sign_tests(testsCorrected)
            
            # find position of observed value in respect to random control
            empiricalPvalue, numberAbove = self.empirical_pvalue(listRandomSignificantsNoWarning, countSignificantNoWarning)
            
            # if empirical value is significant  
            sign = "0"
            if float(empiricalPvalue) < EnrichmentAnalysisStrategy.SIGN_VALUE_AGAINST_RANDOM_CONTROL:
                sign = "1"
                
            outHandlerStats.write("%s\t%i\t%.2f\t%i\t%.1e\t%s\n" % (rnaID, countSignificantNoWarning, avgSignRandomNoWarning, numberAbove, empiricalPvalue, sign))


        # Logger.get_instance().info( "EnrichmentAnalysisStrategy.enrichement_analysis: Tests performed: %s " % str( performedTests ) )
        # Logger.get_instance().info( "EnrichmentAnalysisStrategy.enrichement_analysis: Tests skipped: %s " % str( skippedTests ) )
   
        outHandler.close()
        outHandlerStats.close()


    # #
    # Randomize protein-RNA interactions. 
    # Number of interactions per protein is kept, and number of interactions per lncRNa is kept, only the pairwise combinations are randomised.
    # Not controlling that all interactions are different from real (some may be same as real), so that it works for saturation conditions.
    def randomize_interactions(self, interactions):

        # first filter out interactions with proteins that have no annotation
        filtInteractions = []
        for inter in interactions:
            protID = str(inter.proteinID)
            # only store info of proteins that have annotation information
            if protID not in self.protAnnotDict:
                continue            
            filtInteractions.append( inter)
                                
        # real list of RNA interactions (with duplicates) to maintain same number of interactions per RNA
        listRNAInteractions = [ str(inter.transcriptID) for inter in filtInteractions] 

        #=================================================================== 
        # Randomise protein-RNA interactions
        #=================================================================== 

        # For each RNA, store all proteins it interacts with and their annotations
        randomInteractions = {}  # Key -> randomisation ID, value -> dict; key -> transcript ID, value -> dict; key -> pathway ID, value -> list of prot IDs (after filtering)

        for i in xrange(0, self.numberRandomizations):

            # initialise randomisation
            if i not in randomInteractions:
                randomInteractions[ i] = {}
                
            randomPairs = set() # stores already used interactions

            # shuffle order of appearance of RNAs in each interaction
            randListRNAInteractions = random.sample( listRNAInteractions, len( listRNAInteractions))

            for j in xrange(0, len(filtInteractions)):
                inter = filtInteractions[ j]
                                
                protID = str(inter.proteinID)

                # pick randomised RNA
                txID = randListRNAInteractions[ j]

                ## force all interactions to be different (cannot have duplicated interactions)
                # search other RNA interactor, until one that is not yet used for this protein
                tag = txID+"|"+protID                 
                k = 0 
                while tag in randomPairs:
 
                    k += 1
 
                    try: # it can happen that a final unique interaction is not found, but in large datasets its probability is rare
                        txID = randListRNAInteractions[ j + k]
                    except IndexError:
                        Logger.get_instance().warning("EnrichmentAnalysisStrategy.randomize_interactions : randomized protein-RNA interaction is repeated. %s" % (tag))           
                        break # keep last txID (duplicated interaction)
                    tag = txID+"|"+protID
 
                # replace random list, if possible
                try:
                    randListRNAInteractions[ j + k] = randListRNAInteractions[ j]
                except IndexError:
                    pass
 
                randomPairs.add( tag)
                
                if txID not in randomInteractions[ i]:
                    randomInteractions[ i][ txID] = {}

                for annot in self.protAnnotDict[ protID]:
                    if annot not in randomInteractions[ i][ txID]:
                        randomInteractions[ i][ txID][ annot] = []
                    randomInteractions[ i][ txID][ annot].append(protID)
            
        return randomInteractions
        

    # #
    # Function to run hypergeometric test of RNA against list of annotations
    def run_rna_vs_annotations(self, rna_id, annotation_dict, total_rna_interactions):
        
        pvalues = numpy.empty(len(annotation_dict)) * numpy.nan  # stores list of pvalues to be corrected                  
        tests = numpy.empty(len(annotation_dict), object)  # stores output from tests
                                
        counter = 0
        # for each annotation with at least one interacting partner
        for annotID in annotation_dict:

            # keep this as list so that it works even if reshuffling and having replicate proteins in same annotation.
            protList = []            
            for prot in annotation_dict[ annotID]:
                if prot in total_rna_interactions:
                    protList.append(prot)
            
            # number of proteins with annotation that have interaction predictions
            
            possibleProtList = self.annotWithInteractionDict[ annotID]
            
            assert (len(protList) <= len(possibleProtList)), "RNA cannot interact with more proteins of annotation than existing in annotation"

            # tag whether test does not pass the minimum number of annotations or interactions
            skipTest = 0

            # Skip test/output if minimum number of items is not respected                   
            if len(possibleProtList) < self.minimumProteinAnnotation or len(protList) < self.minimumProteinInteraction :
                skipTest = 1

            #===================================================================   
            # Expression warning 
            #===================================================================          

            # Skip test/output if group expression filter is on and filter is not passed
            if self.expressionWarning != OptionConstants.DEFAULT_EXPRESSION_WARNING:
                tissueCounts = {}  # key -> tissue name, val -> count of proteins present in tissue
                for tiss in self.expressionDict[ rna_id]: 
                    if tiss not in tissueCounts:
                        tissueCounts[ tiss] = 0

                # for each interacting protein, check which tissues they are expressed
                for prot in protList:
                    if prot in self.protTissueExpressions:
                        for tiss in self.expressionDict[ rna_id]: 
                            if tiss in self.protTissueExpressions[ prot]:
                                # i.e. if both RNA and protein present in this tissue
                                tissueCounts[ tiss] += 1
                    else:
                        Logger.get_instance().warning("run_rna_vs_annotations : protein without expression data: %s" % (prot))

                maxTissueOverlap = 0  # store maximum number of proteins present in tissue
                for tiss in tissueCounts:
                    if tissueCounts[ tiss] > maxTissueOverlap:
                        maxTissueOverlap = tissueCounts[ tiss]

                assert maxTissueOverlap <= len(protList), "Maximum number of proteins with interactions in annotation in the same tissue must always be less than number of interacting proteins in annotation"

                # calculate maximum proportion of proteins in same tissue
                if len(protList) != 0:
                    overlapProportion = float(maxTissueOverlap) / float(len(protList))
                else:
                    # division by zero
                    overlapProportion = 0

                # if the proportion in same tissue is not sufficient compared to user input proportion
                if overlapProportion < self.expressionWarning:
                    skipTest = 1


            #===================================================================   
            # Values for hypergeometric test
            #===================================================================          

            # Hypergeometric test parameters
            # note: test for each RNA - annotation pair. 
            x = len(protList)  # white balls drawn ( proteins with current annotation with positive interactions)
            m = len(possibleProtList)  # total white balls ( proteins with the current annotation (that have interaction predictions))
            n = len(self.backgroundProteins) - m  # total black balls ( all the proteins with at least one positive interaction and at least one annotation, but not in current annotation)
            k = len(total_rna_interactions)  # total number of draws ( proteins with positive interactions, regardless of current annotation)

            assert (m + n >= k), "number of draws should be less than total number of balls"

            # code speed up, don't need to recalculate hypergeom test if already done previously
            tag = "%s,%s,%s,%s" % (x, m, n, k)
            if tag in self.testContainer:
                hyperResult = self.testContainer[ tag]
            else:
                # perform the actual test
                hyperResult = self.hypergeometric_test(x, m, n, k)  # this is the slow part of the code
                self.testContainer[ tag] = hyperResult

            self.countTotalTests += 1

            # store pvalue to be corrected                
            pvalues[ counter] = hyperResult               

            percentInteracting = "%.1f%%" % (k * 100.0 / len(self.backgroundProteins))

            text = "%s\t%s\t%s\t%s\t%s\t%s\t%.1e" % (rna_id, annotID, x, m, percentInteracting, skipTest, hyperResult) 
            tests[ counter] = text.split("\t") 
 
            counter += 1
 
        #===================================================================   
        # For each test performed with this RNA
        #===================================================================          
        # calculate corrected p values
                
        nTests = len(annotation_dict)
                
        testsCorrected = self.correct_pvalues(nTests, pvalues, tests)

        return testsCorrected


    # #
    # Function to correct pvalues, and add correction and significance to provided 'tests' list
    def correct_pvalues(self, nTests, pvalues, tests):
        
        testsCorrected = numpy.empty(nTests, object)  # stores data plus correction
         
        correctedPvalues = self.multiple_test_correction(pvalues)
        
        assert (len(pvalues) == len(correctedPvalues))

        # Correction of test selection, by ranking original pvalues, apply correction, 
        # select all the original p-values below the first where corrected pvalue is below our threshold

        # sort pvalues and keep index        
        sortedPvalues = sorted(pvalues, reverse=True)
        sortedPvaluesIdx = sorted(xrange(len(pvalues)), reverse=True, key=lambda k: pvalues[ k])

        # correspondece of pvalues -> corrected pvalues of the sorted pvalues array
        correctedPvaluesOfSorted = [ correctedPvalues[ idx]  for idx in sortedPvaluesIdx]

        # stores list of indexes deemed significant with our own selection of rejected values from 
        significantIndexes = []

        # find the first corrected pvalue below our threshold, and pick all the original pvalues on the sorted array as accepted            
        for i in xrange(len(sortedPvalues)):
            if correctedPvaluesOfSorted[ i] < EnrichmentAnalysisStrategy.SIGN_VALUE_TEST:
                significantIndexes = sortedPvaluesIdx[ i: ]
                break

#         # to observe that some corrected pvalues will be higher than threshold, yet their original pvalue is below the first sorted pvalue to be accepted
#         for idx in significantIndexes:
#              
#             if correctedPvalues[ idx] > EnrichmentAnalysisStrategy.SIGN_VALUE_TEST:
#                 print "SEE", correctedPvalues[idx], pvalues[idx]

        for i in xrange(0, len(tests)):
            l = tests[i][:]
   
            # add corrected pvalue to existing list
            corr = "%.1e" % correctedPvalues[i]
            l.append(corr)
   
            # significative result tag
            sign = "0"
            if i in significantIndexes:
                sign = "1"
   
            # add sign tag to existing list
            l.append(sign)
   
            # append current list to list of RNA vs annotation
            testsCorrected[ i] = l

        return testsCorrected

    # #
    # Run multipletest correction and return pvalues
    def multiple_test_correction(self, pvalues, meth="fdr_bh"):
        return multipletests(pvalues, method=meth)[1]


    # #
    # Run hypergeometric test using scipy. Based on R phyper rationel.
    def hypergeometric_test(self, x, m, n, k):

        # Documentation from R phyper function.
        # x, q vector of quantiles representing the number of white balls drawn
        # without replacement from an urn which contains both black and white
        # balls.
        # 
        # m the number of white balls in the urn.
        # 
        # n the number of black balls in the urn.
        # 
        # k the number of balls drawn from the urn.

        # Notes on behaviour of the test:
        # - If all white balls are found (x == m), pvalue will be 0.0, regardless of other parameters
        # - If there are no white balls (m), pvalue will be 0.0, regardless of other parameters
        # - If there are no balls drawn (k), pvalue will be 0.0, regardless of other parameters
        # - If there are not white balls withdrawn (x), but number of m and k is small , p-val may still be significant e.g. x = 0, m = 1, n = 36, k = 1. pval = 0.027

#         cmd = "Rscript /home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/Rscripts/hypergeom_test.R %s %s %s %s" % (x,m,n,k)
#         os.system(cmd)
        
        if k > (m + n):
            raise RainetException("EnrichmentAnalysisStrategy.hypergeometric_test: Number of draws is higher than total balls.")                   
        if m == 0:
            raise RainetException("EnrichmentAnalysisStrategy.hypergeometric_test: Number of white balls is zero.")
        if k == 0:
            raise RainetException("EnrichmentAnalysisStrategy.hypergeometric_test: Number of draws is zero.")
        if n < 0:
            raise RainetException("EnrichmentAnalysisStrategy.hypergeometric_test: Number of black balls is negative. Possible issue with protein background.")
            
 
        # stats.hypergeome.sf gives the same result as R phyper
        if self.lowerTail:
            testResult = stats.hypergeom.sf(x, m + n, m, k)
        else:
            # equivalent to R lower.tail = FALSE
            testResult = 1 - stats.hypergeom.sf(x, m + n, m, k)
            

        # if x == 0 and testResult <= EnrichmentAnalysisStrategy.SIGN_VALUE:
        #    Logger.get_instance().warning( "EnrichmentAnalysisStrategy.hypergeometric_test: Number of white balls withdrawn is zero but result is significant. x=%s, m=%s, n=%s, k=%s" % (x,m,n,k)  )

        # print ("x: %i\tm: %i\tn: %i\tk: %i" % (x,m,n,k) )
        # print ("Hypergeometric_test p-value:\t%.3f\n" % testResult )
            
        return testResult


    # # 
    # Count number of significant tests among tests performed
    def count_sign_tests(self, tests_corrected):

        signColumn = EnrichmentAnalysisStrategy.SIGN_COLUMN
        warningColumn = EnrichmentAnalysisStrategy.WARNING_COLUMN

        countSignificant = 0  # number of significant tests regardless of warning
        countSignificantNoWarning = 0  # number of significant tests without a warning
        for test in tests_corrected:

            # if test is significant and does not have warning flag
            if test[ signColumn] == "1":
                countSignificant += 1
            elif test[ signColumn] == "0":
                pass
            else:
                raise RainetException("EnrichmentAnalysisStrategy.count_sign_tests: Significance value is not boolean.")

            if test[ signColumn] == "1" and test[ warningColumn] == "0":
                countSignificantNoWarning += 1

        return countSignificant, countSignificantNoWarning


    # #
    # Calculate proportion of random tests that are below observed value
    # Conservative approach: only counts observed below if < random value (not <=)
    def empirical_pvalue(self, list_random_significant, observed_significant, aboveTail=1):
        
        listRandomSignificant = numpy.sort(list_random_significant)            

        below = 0  # count how many times observed value is below random values
        above = 0  # count how many times observed value is above random values
        for val in listRandomSignificant:
            if observed_significant < val:  # conservative
                below += 1
            if observed_significant > val:
                above += 1

        if aboveTail:
            pval = 1.0 - ((above + 1) / float(len(listRandomSignificant) + 1))
            return pval, above
        else:
            pval = 1.0 - ((below + 1) / float(len(listRandomSignificant) + 1))  
            return pval, below


    # #
    # Randomize values in a dictionary while keeping the structure of the dictionary
    def randomize_annotation(self, annotDict, listOfProteins):
                
        # shuffle list of proteins (use of sample with maximum number of sample size, same as shuffle)
        randomizedListOfProteins = random.sample(listOfProteins, len(listOfProteins))
        
#         randomizedListOfProteins = listOfProteins[:]
#         random.shuffle(randomizedListOfProteins)
        
        assert len(randomizedListOfProteins) == len(listOfProteins)

        assert len(annotDict) > 1, "randomization only works if more than one annotation"

        # container of randomized annotation dict
        randomAnnotDict = {}  # key -> annot, val -> list of proteins

        counter = 0
        
        # for each annotation of original annotation dict, add proteins randomly but with same annotation topology
        for annot in sorted(annotDict):  # sorted is important
            
            if annot not in randomAnnotDict:
                randomAnnotDict[ annot] = []
            
            nItems = len(annotDict[ annot])

            # for each protein in annotation            
            for i in xrange(0, nItems):
                randomAnnotDict[ annot].append(randomizedListOfProteins.pop())

#                 currentProt = randomizedListOfProteins[ counter]
#                  
#                 # Note: if using NetworkModules, there are overlapping annotations.
#                 # I.e. same protein present in different annotations
#                 # avoid replicates of same protein falling into same annotation
#                 # if protein already in bag. Use of while in case this happens consecutively  
#                 if currentProt in randomAnnotDict[ annot]:
# 
#                     if len( randomAnnotDict) == 1:
#                         # Approach 1: in case there is only one annotation initialised in dictionary, swap with next protein
# 
#                         swaps = 1
#                         # if protein already in bag. Use of while in case this happens consecutively
#                         while currentProt in randomAnnotDict[ annot]:
#            
#                             otherProt = randomizedListOfProteins[counter + swaps]
# 
#                             if otherProt != currentProt:
#                                 # swap current with next
#                                 toBeSwapped = randomizedListOfProteins[counter ]
#                                 randomizedListOfProteins[ counter] = otherProt
#                                 randomizedListOfProteins[ counter + swaps] = toBeSwapped
#                
#                                 currentProt = randomizedListOfProteins[ counter]
#                                 swaps+=1
# 
#                     else:
#                         # Approach 2: swap with protein of a previously filled annotation
#                      
#                         succ = 0
#                         tries = 0
#       
#                         # steal protein from another bag   
#                         while succ == 0:
#                             for annot2 in randomAnnotDict:                           
#                                  
#                                 if annot == annot2:
#                                     continue
#      
#                                 otherProt = randomAnnotDict[ annot2][ tries]
#                                  
#                                 # if the new protein is different than current, and the current is not already in bag from where protein is being stolen
#                                 if otherProt != currentProt and currentProt not in randomAnnotDict[ annot2]:
#                                     # swap current with other
#                                     toBeSwapped = currentProt
#                                     randomAnnotDict[ annot2][ 0] = toBeSwapped
#                                     currentProt = otherProt
#                                     # successful swap
#                                     succ = 1
#                                     break
#                                 else: # go to next bag
#                                     continue
#      
#                             tries += 1
#      
#                         assert succ == 1
# 
#                       
#                 randomAnnotDict[ annot].append( currentProt)
#                   
#                 counter += 1


        assert len(randomAnnotDict) == len(annotDict)
            
        return randomAnnotDict


    # #
    # Randomize values in a dictionary while keeping the structure of the dictionary
    # Approach: swap the identify of the protein (e.g. Protein1 becomes Protein2, Protein5 becomes Protein1 etc)
    def randomize_proteins(self, annotDict, listOfProteins):

        # remove redundancy in list           
        listOfProteins = list(set(listOfProteins))

        # shuffle list of proteins (use of sample with maximum number of sample size, same as shuffle)
        randomizedListOfProteins = random.sample(listOfProteins, len(listOfProteins))

        # assert sorted( listOfProteins) == sorted( randomizedListOfProteins)

        # assign a number to each protein
        proteinShapeshifter = {}  # key -> old ID of protein, val -> new ID of protein
        for i in xrange(len(listOfProteins)):
            oldProtID = listOfProteins[ i]
            if oldProtID not in proteinShapeshifter:
                proteinShapeshifter[ oldProtID] = ""
            else:
                raise RainetException("EnrichmentAnalysisStrategy.randomize_proteins: randomization issue: original ID is not unique", listOfProteins[ i])
            
            # assign new identity to protein
            proteinShapeshifter[ oldProtID] = randomizedListOfProteins[ i]

        assert len(proteinShapeshifter) == len(listOfProteins)

        # container of randomized annotation dict
        randomAnnotDict = {}  # key -> annot, val -> list of proteins

        # for each annotation of original annotation dict, change the protein IDs to the swapped ones
        for annot in sorted(annotDict):  # sorted is important
             
            if annot not in randomAnnotDict:
                randomAnnotDict[ annot] = []
 
            for prot in annotDict[ annot]:
                newProt = proteinShapeshifter[ prot]
                randomAnnotDict[ annot].append(newProt)
 
        assert len(randomAnnotDict) == len(annotDict)

        return randomAnnotDict


    # #
    # Retrieve expression per tissue for each RNA and protein.
    # needs access to list of proteins with annotation, therefore it needs to be run after self.annotation_report
    def retrieve_expression(self):

        self.expressionDict = {}  # key -> transcript ID, value -> set of tissues passing expression cutoff

        self.protTissueExpressions = {}  # key -> protein ID, value -> set of tissues passing expression cutoff

        #===================================================================    
        # Map mRNA to protein ID
        #===================================================================    
        # Search mRNA that produces interacting protein
        mRNAMap = self.sql_session.query(MRNA.transcriptID, MRNA.proteinID).all()
        mRNADict = {}  # key -> protein ID, val -> list of mRNAs encoding protein
        for items in mRNAMap:
            txID, protID = str(items[0]), str(items[1])
            if protID != "None":
                if protID not in mRNADict:
                    mRNADict[ protID] = []
                mRNADict[ protID].append(txID)
        
        Logger.get_instance().info("retrieve_expression : initialised mRNA-protein data. %s proteins with mRNAs." % len(mRNADict))    


        #===================================================================    
        # Map expression per tissue to transcript ID   
        #===================================================================    
        # Load entire expression table        
        expressionMap = self.sql_session.query(RNATissueExpression.transcriptID, RNATissueExpression.expressionValue, RNATissueExpression.tissueName).all()

        # store tissues where present for each RNA
        for items in expressionMap:
            txID, expr, tissName = str(items[0]), float(items[1]), str(items[2])
            if txID not in self.expressionDict:
                self.expressionDict[ txID] = set()
            
            if expr >= self.minimumExpression:
                self.expressionDict[ txID].add(tissName)
        
        Logger.get_instance().info("retrieve_expression : loaded expression data. %s total RNAs with expression data loaded." % len(self.expressionDict))    
        
        #===================================================================    
        # Store all protein-related expression values into memory
        #===================================================================    
        self.protTissueExpressions = {}  # key -> prot ID, val -> dict. key -> tissues where present
        for protID in self.proteinWithAnnotationWithInteraction:
            # skip protein with no mRNAs in database
            if protID not in mRNADict:
                continue
        
            if protID not in self.protTissueExpressions:
                self.protTissueExpressions[ protID] = set()
        
            # Search mRNA that produces interacting protein
            mRNAs = mRNADict[ protID]
        
            # Get Protein expression for all tissues
            # there can be several mRNAs for the same protein ID, here we use them all to have set of tissues where present
            # i.e. we only require that at least one of the mRNAs producing the protein is present in a tissue
            for mRNAID in mRNAs:
            
                # skip transcripts with no expression      
                if mRNAID not in self.expressionDict:
                    continue
                                                  
                # Get the Protein expression, using its mRNA
                for tissue in self.expressionDict[ mRNAID]:                         
                    self.protTissueExpressions[ protID].add(tissue)
        
        Logger.get_instance().info("retrieve_expression : initialised expression data. %s proteins with expression data loaded." % len(self.protTissueExpressions))    


    # #
    # Run Rscript to produce Sweave file and consequent pdf report, using the data written by this script
    def write_report(self):
        pass

