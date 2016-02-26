
import os
#import pandas as pd
from sqlalchemy import or_,and_,distinct

from fr.tagc.rainet.core.execution.ExecutionStrategy import ExecutionStrategy
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.data.DataManager import DataManager

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

from fr.tagc.rainet.core.data.TableStatus import TableStatus
from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util import Constants
from fr.tagc.rainet.core.data.RNATissueExpression import RNATissueExpression


# #
# This class define the Strategy enabling user to perform DB query interactively
class AnalysisStrategy(ExecutionStrategy):


    #===============================================================================
    # Analysis strategy Constants
    #===============================================================================

    # Data Manager Keywords

    RNA_ALL_KW = "allRNAs"
    PROT_ALL_KW = "allProteins"

    RNA_FILTER_KW = "selectedRNAs"
    PROT_FILTER_KW = "selectedProteins"
    PRI_FILTER_KW = "selectedInteractions"

    # Report files
    PARAMETERS_LOG = "parameters.log"
    REPORT_RNA_NUMBERS = "rna_numbers.tsv"
    REPORT_RNA_EXPRESSION = "rna_expression.tsv"
    REPORT_RNA_EXPRESSION_DATA_PRESENCE = "rna_expression_data_presence.tsv"
    REPORT_INTERACTION_NUMBERS = "interaction_numbers.tsv"

    # #
    # The Strategy execution method
    def execute(self):

        #===================================================================
        # Getting input arguments        
        #===================================================================
        
        self.DBPath = OptionManager.get_instance().get_option(OptionConstants.OPTION_DB_NAME)
        self.species = OptionManager.get_instance().get_option(OptionConstants.OPTION_SPECIES)
        self.outputFolder = OptionManager.get_instance().get_option(OptionConstants.OPTION_OUTPUT_FOLDER)
        self.minimumInteractionScore =  OptionManager.get_instance().get_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE)
        self.transcriptBiotype = OptionManager.get_instance().get_option(OptionConstants.OPTION_TRANSCRIPT_BIOTYPE)
        self.lncRNABiotypes = OptionManager.get_instance().get_option(OptionConstants.OPTION_LNCRNA_BIOTYPES)
        self.gencode = OptionManager.get_instance().get_option(OptionConstants.OPTION_GENCODE)
        self.expressionValueCutoff = OptionManager.get_instance().get_option(OptionConstants.OPTION_EXPRESSION_VALUE_CUTOFF)

        self.arguments = {OptionConstants.OPTION_DB_NAME : self.DBPath,
                          OptionConstants.OPTION_SPECIES : self.species,
                          OptionConstants.OPTION_OUTPUT_FOLDER : self.outputFolder,
                          OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE : self.minimumInteractionScore,
                          OptionConstants.OPTION_TRANSCRIPT_BIOTYPE : self.transcriptBiotype,
                          OptionConstants.OPTION_LNCRNA_BIOTYPES : self.lncRNABiotypes,
                          OptionConstants.OPTION_GENCODE : self.gencode,
                          OptionConstants.OPTION_EXPRESSION_VALUE_CUTOFF : self.expressionValueCutoff
                        }

        
        # Check if output folder path exists, create it if not
        if self.outputFolder != "" and len(self.outputFolder) > 0:
            FileUtils.initialise_output_folders(self.outputFolder)
            self.outputFolderReport = self.outputFolder+"/"+Constants.REPORT_FOLDER
        else:
            raise RainetException( "AnalysisStrategy.execute: Provided output folder is empty.")
            
        
        # Check if minimumInteractionScore is float or OFF
        if self.minimumInteractionScore != OptionConstants.DEFAULT_INTERACTION_SCORE:
            try:
                float(self.minimumInteractionScore)
            except TypeError:
                raise RainetException( "AnalysisStrategy.execute: Provided minimum interaction score is not a float.")

        # Check if transcript biotype is correct
        if self.transcriptBiotype not in OptionConstants.TRANSCRIPT_BIOTYPES:
            raise RainetException( "AnalysisStrategy.execute: Provided transcript biotype is not allowed: "+str(self.transcriptBiotype))

        # Process and check provided list of lncRNA subtypes
        if self.lncRNABiotypes != OptionConstants.DEFAULT_LNCRNA_BIOTYPES:
            bioTypeStr = ""
            for subtype in self.lncRNABiotypes.split(","):
                if subtype not in OptionConstants.LNCRNA_BIOTYPES:
                    raise RainetException( "AnalysisStrategy.execute: Provided lncRNA biotype is not allowed: "+str(subtype))              
                bioTypeStr+='"'+subtype+'",'
            self.lncRNABiotypes = bioTypeStr[:-1] # taking out last comma

        # Check if gencode argument is correct
        try:
            self.gencode = int(self.gencode)
        except TypeError:
            raise RainetException( "AnalysisStrategy.execute: Provided GencodeBasicOnly argument is not numeric.")
        if self.gencode != 1 and self.gencode != 0:
            raise RainetException( "AnalysisStrategy.execute: Provided GencodeBasicOnly argument must be either 0 or 1.")


        #===================================================================
        # Initialisation
        #===================================================================

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.DBPath)
        self.sql_session = SQLManager.get_instance().get_session()
        self.db_engine = SQLManager.get_instance().get_engine()
                        
        self.analysis()
        
    # #
    # Run analysis-related functions in order
    def analysis(self):

        Logger.get_instance().info( "AnalysisStrategy.analysis: Starting..." )
        
        #===================================================================
        # Filter datasets
        #===================================================================

        self.filter_RNA()
        
        self.filter_protein()
        
        self.filter_PRI()

        #===================================================================
        # Perform analysis
        #===================================================================

        self.after_filter_report()
# 
#         self.enrichement_analysis()


    # #
    # Filter RNA models
    #
    # Stores final list of RNAs on DataManager 
    def filter_RNA(self):

        Logger.get_instance().info("AnalysisStrategy.filter_RNA..")

        #===================================================================
        # Get all RNA objects
        #===================================================================
        allRNAs = self.sql_session.query(RNA).all()
        DataManager.get_instance().store_data(AnalysisStrategy.RNA_ALL_KW, allRNAs)

        #===================================================================              
        # Filter transcripts based on biotype or Filter LncRNAs based on wanted subtypes (if lncRNA biotype chosen)
        #===================================================================
        if self.transcriptBiotype == OptionConstants.BIOTYPE_LNCRNA and self.lncRNABiotypes != OptionConstants.DEFAULT_LNCRNA_BIOTYPES:       
            queryText = "query(" + self.transcriptBiotype + ".transcriptID).filter(" + self.transcriptBiotype + ".transcriptBiotype.in_([" + self.lncRNABiotypes + "])).all() "
        else:
            queryText = "query(" + self.transcriptBiotype + ".transcriptID).all()"

        filterRNA1 = eval('self.sql_session.' +  queryText)
        filterRNA1 = {str(item.transcriptID) for item in filterRNA1}

        #===================================================================
        # Filter transcripts based on gencode_basic presence
        #===================================================================
        if self.gencode == 1:
            queryText = "query(RNA.transcriptID).filter(RNA.transcriptGencodeBasic == 1).all()"
        else:
            queryText = "query(RNA.transcriptID).all()"

        filterRNA2 = eval('self.sql_session.' +  queryText)
        filterRNA2 = {str(item.transcriptID) for item in filterRNA2}

        #===================================================================                               
        # Get intersection of the various filterings
        # Store into data manager
        #===================================================================
        selectedRNAs = []
        for rna in allRNAs:
            if rna.transcriptID in filterRNA1 and rna.transcriptID in filterRNA2:
                selectedRNAs.append(rna)

        DataManager.get_instance().store_data(AnalysisStrategy.RNA_FILTER_KW, selectedRNAs)


        # TODO: Select one transcript isoform per gene


    # #
    # Filter protein models    
    #
    # Stores final list of proteins on DataManager 
    def filter_protein(self):

        Logger.get_instance().info("AnalysisStrategy.filter_protein..")
        
        #===================================================================
        # Get all Protein objects
        #===================================================================        

        allProts = "query( Protein ).all()"
        DataManager.get_instance().store_data(AnalysisStrategy.PROT_ALL_KW, allProts)

        # TODO: potential filterings on the protein level
        
        selectedProts = allProts # to be changed

        DataManager.get_instance().perform_query(AnalysisStrategy.PROT_FILTER_KW, selectedProts) 


    # #
    # Filter protein-RNA interactions
    def filter_PRI(self):

        Logger.get_instance().info("AnalysisStrategy.filter_PRI..")

        #===================================================================
        # Retrieve selected RNAs and proteins
        #===================================================================                       
        selectedRNAs = {str(item.transcriptID) for item in DataManager.get_instance().get_data(AnalysisStrategy.RNA_FILTER_KW) }
        selectedProteins = {str(item.uniprotAC) for item in DataManager.get_instance().get_data(AnalysisStrategy.PROT_FILTER_KW) } 

        #===================================================================         
        # Filter interactions based on minimumInteractionScore
        #===================================================================                 
        if self.minimumInteractionScore != OptionConstants.DEFAULT_INTERACTION_SCORE:
            queryText = "query( ProteinRNAInteractionCatRAPID ).filter(ProteinRNAInteractionCatRAPID.interactionScore >= "+str(self.minimumInteractionScore)+").all()"    
        else:
            # Running this on whole database may crash computer
            nItems = self.sql_session.query( ProteinRNAInteractionCatRAPID ).count()
            if nItems > 1000000: # dr: improve this
                raise RainetException( "AnalysisStrategy.filter_PRI: intended query may use prohibitive amounts of system memory. Quitting..")
            else:
                queryText = "query( ProteinRNAInteractionCatRAPID ).all()"
 
        interactions = eval('self.sql_session.' +  queryText)

        #===================================================================          
        # Filter for interactions between selected RNAs and proteins
        #===================================================================         
        selectedInteractions = []
        for inter in interactions:
            if inter.transcriptID in selectedRNAs and inter.proteinID in selectedProteins:
                selectedInteractions.append(inter)
 
        DataManager.get_instance().store_data(AnalysisStrategy.PRI_FILTER_KW, selectedInteractions)

        #===================================================================    
        # Filter interactions based on RNA and Protein (mRNA) expression
        #=================================================================== 
        # TODO: first we need to know which expression data we have


    # #
    # Retrieve statistics before and after the filtering steps.
    # Produce output files that will be used for a pdf report
    def after_filter_report(self):

        #===================================================================    
        # Write log of parameters used
        #=================================================================== 

        outHandler = FileUtils.open_text_w( self.outputFolder + "/" + AnalysisStrategy.PARAMETERS_LOG )
        
        Logger.get_instance().info( "\nARGUMENTS USED:" )

        for argName,argValue in self.arguments.iteritems():
            Logger.get_instance().info( "%s: %s" % ( argName, argValue) ) 
            outHandler.write( "%s: %s\n" % ( argName, argValue) )
        outHandler.close()

        #===================================================================    
        # RNA numbers report
        #
        # 'Before filtering' = before any RNA, protein or interactions filtering
        # 'After filtering' = after RNA filter, but before interactions filter
        #=================================================================== 

        #
        # File with number of RNAs types and lncRNA subtypes, before and after filter
        #

        outHandler = FileUtils.open_text_w( self.outputFolderReport + "/" + AnalysisStrategy.REPORT_RNA_NUMBERS )
        
        # Write header
        outHandler.write( "Data\t" +
                          "Gene\t" +
                           "\t".join( OptionConstants.TRANSCRIPT_BIOTYPES) + "\t" +
                           "\t".join( OptionConstants.LNCRNA_BIOTYPES) + "\n")

        # Get filtered and unfiltered RNAs
        filteredRNAs = DataManager.get_instance().get_data( AnalysisStrategy.RNA_FILTER_KW)
        allRNAs = DataManager.get_instance().get_data( AnalysisStrategy.RNA_ALL_KW)

        # Get their types
        allRNATypes = [ rna.type for rna in allRNAs]
        filteredRNATypes = [ rna.type for rna in filteredRNAs]

        allLncRNATypes = [ rna.transcriptBiotype for rna in allRNAs if rna.type == OptionConstants.BIOTYPE_LNCRNA]
        filteredLncRNATypes = [ rna.transcriptBiotype for rna in filteredRNAs if rna.type == OptionConstants.BIOTYPE_LNCRNA]

        allGenes = { rna.geneID for rna in allRNAs}
        filteredGenes = { rna.geneID for rna in filteredRNAs}

        # # Report numbers before and after filtering        
        beforeFilterText = "Before_RNA_filter"
        afterFilterText = "After_RNA_filter"

        # Number of unique gene IDs
        beforeFilterText+= "\t%i" % len(allGenes)
        afterFilterText+= "\t%i" % len(filteredGenes)      
        
        # RNA biotypes
        for biotype in OptionConstants.TRANSCRIPT_BIOTYPES:
            if biotype != OptionConstants.DEFAULT_BIOTYPE:
                beforeFilterText+= "\t%i" % allRNATypes.count( biotype)
                afterFilterText+= "\t%i" % filteredRNATypes.count( biotype)
            else:
                beforeFilterText+= "\t%i" % len( allRNAs)
                afterFilterText+= "\t%i" % len( filteredRNAs)
        
        # LncRNA biotypes
        for lncBiotype in OptionConstants.LNCRNA_BIOTYPES:
            beforeFilterText+=  "\t%i" % allLncRNATypes.count( lncBiotype)
            afterFilterText+=  "\t%i" % filteredLncRNATypes.count( lncBiotype)

        outHandler.write( beforeFilterText+"\n"+afterFilterText+"\n")
        outHandler.close()

        #===================================================================    
        # RNA expression report
        #=================================================================== 

        #
        # File with average expression (among tissues) for each transcript, discrimination of RNA types and lncRNA subtypes
        #

        outHandler = FileUtils.open_text_w( self.outputFolderReport + "/" + AnalysisStrategy.REPORT_RNA_EXPRESSION )
        
        # Write header
        outHandler.write("transcriptID\ttype\ttranscriptBiotype\tmeanExpression\n") 

        expressionDataCounts = {} # stores counts of RNAs with or without expression data per RNA subtype
        for rna in filteredRNAs:
            # Get expression values on the several tissues for filtered transcript
            queryResult = self.sql_session.query( RNATissueExpression.expressionValue).filter( RNATissueExpression.transcriptID == rna.transcriptID).all()

            # Initialise dict to store existence/absence of data
            if rna.transcriptBiotype not in expressionDataCounts:
                expressionDataCounts[ rna.transcriptBiotype] = {}
                expressionDataCounts[ rna.transcriptBiotype]["with"] = 0
                expressionDataCounts[ rna.transcriptBiotype]["without"] = 0

            # if there is expression data for this transcript
            if queryResult != None and len(queryResult) > 0:
                expressionDataCounts[ rna.transcriptBiotype]["with"]+= 1

                # Array with the actual expression values
                expressionValues = [ result[0] for result in queryResult]        
                # Write into file the average expression value between tissues
                outHandler.write( "%s\t%s\t%s\t%.2f\n" % (rna.transcriptID, rna.type, rna.transcriptBiotype, sum(expressionValues)/len(expressionValues)) )
            else:
                # If is possible to have transcripts (from RNA table) that are not present in the "RNATissueExpression" table,
                # since we and GTEx are using different Ensembl/GENCODE releases and some transcripts were deprecated or are new.
                expressionDataCounts[ rna.transcriptBiotype]["without"]+= 1

        outHandler.close()

        #
        # File with percentage of transcript with expression data, discrimination of RNA types and lncRNA subtypes
        #

        outHandler = FileUtils.open_text_w( self.outputFolderReport + "/" + AnalysisStrategy.REPORT_RNA_EXPRESSION_DATA_PRESENCE )

        # Write header
        outHandler.write("subtype\ttx_with_expression_data\ttx_without_expression_data\tperc_tx_with_expression_data\n") 

        for subtype in expressionDataCounts:
            withExpression = expressionDataCounts[subtype]["with"]
            withoutExpression = expressionDataCounts[subtype]["without"]
            if withExpression > 0:
                perc = "%.2f%%" % (withExpression*100.0 / (withoutExpression+withExpression) )
            else:
                perc = "NA"
            outHandler.write("%s\t%i\t%i\t%s\n" % ( subtype, withExpression, withoutExpression, perc) )
        
        outHandler.close()


        #===================================================================    
        # Interactions report
        #
        # Note1: interactions filtering is applied after RNA and protein filterings
        #        'Before filtering' = before interactions filtering, but after RNA and protein filters
        #        'After filtering' = after interactions filtering AND after RNA and protein filters
        #        Therefore, 'filtering'-named objects refer to the interaction filter
        #
        # Note2: the protein-RNA interactions are a natural filter, even with default parameters,
        #        not all RNA's / proteins will be interacting (e.g. they lack interaction data)
        #=================================================================== 

        #
        # File with number of interactions and numbers of proteins and RNAs involved, before and after filter
        #

        outHandler = FileUtils.open_text_w( self.outputFolderReport + "/" + AnalysisStrategy.REPORT_INTERACTION_NUMBERS )

        # Write header
        outHandler.write( "Data\t" +
                          "Total_interactions\t" +
                          "Interacting_Proteins\t" +
                          "Interacting_RNAs\t" +
                          "Interacting_LncRNAs\t" +                          
                           "\t".join( OptionConstants.LNCRNA_BIOTYPES) + "\n")

        # Get filtered interactions
        filteredInteractions = DataManager.get_instance().get_data(AnalysisStrategy.PRI_FILTER_KW)

        # Get count of filtered and unfiltered interactions
        # Note: opposed to filtered interactions, it may not be possible to retrieve all interaction objects due to computational constraints
        filteredInteractionsCount = len(filteredInteractions)
        allInteractionsCount = self.sql_session.query( ProteinRNAInteractionCatRAPID ).count()
        
        # Lists of IDs for the filtered interactions
        setRNAs = {inter.transcriptID for inter in filteredInteractions}
        setProteins = {inter.proteinID for inter in filteredInteractions}
        
        # Counts for proteins and RNAs with interactions, before and after filtering
        allDistinctTxCount = self.sql_session.query( ProteinRNAInteractionCatRAPID.transcriptID ).distinct().count()
        allDistinctProtCount = self.sql_session.query( ProteinRNAInteractionCatRAPID.proteinID ).distinct().count()
        filteredDistinctTxCount = len(setRNAs)
        filteredDistinctProtCount = len(setProteins)

        # # Report numbers before and after filtering        
        beforeFilterText = "Before_interactions_filter"
        afterFilterText = "After_interactions_filter"

        # Total number of interactions
        beforeFilterText+= "\t%i" % allInteractionsCount
        afterFilterText+= "\t%i" % filteredInteractionsCount   

        # numbers of proteins in interactions
        beforeFilterText+= "\t%i" % allDistinctProtCount
        afterFilterText+= "\t%i" % filteredDistinctProtCount
        
        # numbers of RNAs in interactions
        beforeFilterText+= "\t%i" % allDistinctTxCount
        afterFilterText+= "\t%i" % filteredDistinctTxCount

        # numbers of lncRNAs in interactions, and their subtypes
        filteredRNATypesCounts = {} # key -> transcriptBiotype, value -> number of interacting RNAs with that transcriptBiotype
        filteredInteractingLncRNAs = {} # key -> transcriptID of interacting lncRNA, value -> RNA object

        # Get values for filtered interactions
        for rna in filteredRNAs: # filtered RNAs is list of RNAs after RNA filter (not related to interactions filter)
            # if the RNA is in the set of interacting RNAs
            if rna.transcriptID in setRNAs:
                if rna.transcriptBiotype not in filteredRNATypesCounts:
                    filteredRNATypesCounts[ rna.transcriptBiotype] = 0
                filteredRNATypesCounts[ rna.transcriptBiotype]+= 1

                # if interacting RNA is a lncRNA
                if rna.transcriptBiotype in OptionConstants.LNCRNA_BIOTYPES: 
                    if rna.transcriptID not in filteredInteractingLncRNAs:
                        filteredInteractingLncRNAs[rna.transcriptID] = rna
                    else:
                        RainetException( "AnalysisStrategy.after_filter_report : abnormal duplicate transcript ID: " +  rna.transcriptID)

        # numbers of lncRNAs           
        # query number of lncRNAs that are interacting             
        allInteractingLncRNACount = self.sql_session.query( ProteinRNAInteractionCatRAPID.transcriptID ).filter( \
                                                            and_(ProteinRNAInteractionCatRAPID.transcriptID == RNA.transcriptID,
                                                            RNA.type == OptionConstants.BIOTYPE_LNCRNA ) ).distinct().count()

        assert( allInteractingLncRNACount <= len(allRNAs),
                 "Number of all lncRNAs must be less or equal to total number of all RNAs.")

        assert( len(filteredInteractingLncRNAs) <= len(filteredRNAs),
                 "Number of filtered lncRNAs must be less or equal to total number of filtered RNAs.")

        beforeFilterText+= "\t%i" % allInteractingLncRNACount
        afterFilterText+= "\t%i" % len( filteredInteractingLncRNAs)

        # numbers for each subtype of lncRNAs
        for lncBiotype in OptionConstants.LNCRNA_BIOTYPES:
            # before filter
            # query for all interactions in PRI table for numbers of RNAs with specific biotype
            allRNACount = self.sql_session.query( ProteinRNAInteractionCatRAPID.transcriptID ).filter( \
                                                  and_(ProteinRNAInteractionCatRAPID.transcriptID == RNA.transcriptID, 
                                                  RNA.transcriptBiotype == lncBiotype ) ).distinct().count()

            # Note: a "count" query will give 0 and not None if nothing found            
                 
            # after filter
            if lncBiotype in filteredRNATypesCounts:
                filteredRNACount = filteredRNATypesCounts[ lncBiotype]
            else:
                filteredRNACount = 0

            assert( allRNACount <= allInteractingLncRNACount,
                    "Number of interacting lncRNAs of a subtype must be less or equal to number of interacting lncRNAs.")
            assert( filteredRNACount <= filteredInteractingLncRNAs,
                    "Number of filtered interacting lncRNAs of a subtype must be less or equal to number of filtered interacting lncRNAs.")


            beforeFilterText+= "\t%i" % allRNACount
            afterFilterText+= "\t%i" % filteredRNACount


        outHandler.write( beforeFilterText + "\n" + afterFilterText + "\n")
        outHandler.close()
           
#         # Percentage of interacting RNAs
#         filteredRNAs = len(DataManager.get_instance().get_data(AnalysisStrategy.RNA_FILTER_KW))
#         print ("%% interacting LncRNAs: %.2f" % (len(setRNAs)*100.0/filteredRNAs) )
#         # Percentage of interacting Proteins
#         filteredProteins = len(DataManager.get_instance().get_data(AnalysisStrategy.PROT_FILTER_KW))
#         print ("%% interacting Proteins: %.2f" % (len(setProteins)*100.0/filteredProteins) )


        
    def enrichement_analysis(self):

        pass

#  
#         # Approach: For each interacting RNA, does it target significantly more proteins in same KEGG pathway
#          
#         # 1) for each KEGG pathway, count how many proteins in it
#         keggFrequencies = {}
#         keggPathways = self.sql_session.query( KEGGPathway ).all()
#         for pathway in keggPathways:
#             keggFrequencies[pathway.keggID] = self.sql_session.query( ProteinKEGGAnnotation ).filter( ProteinKEGGAnnotation.keggPathway_id == pathway.keggID ).count()
#         
#         # 2) for each interacting RNA, retrieve kegg pathways of proteins it interacts with
#         interactions = DataManager.get_instance().get_data(DataConstants.PRI_FILTER_KW)


#         # Get protein IDs and transcript IDs present in filtered interactions
#         for inter in selectedInteractions:
#             print (inter.transcriptID,inter.proteinID)
#              
#             prot = self.sql_session.query( Protein ).filter( Protein.uniprotAC == inter.proteinID ).all()
#  
#             rna = self.sql_session.query( RNA ).filter( RNA.transcriptID == inter.transcriptID ).all()

    def hypergeometric_test(self):
        pass

        #scipy.stats.hypergeom.cdf(x, M, n, N, loc=0)


    # #
    # Run Rscript to produce Sweave file and consequent pdf report, using the data written by this script
    def run_statistics(self):
        pass

#         # confirm that required input files are present
#         for filePath in ProcessGTExData.R_REQUIRED_FILES:
#             if not os.path.exists(filePath):
#                 raise RainetException( "run_statistics : Input file is not present: " + filePath )
#                 
#         # launch the analysis
#         command = "cd " + os.path.dirname(ProcessGTExData.SWEAVE_R_SCRIPT) + "; Rscript %s %s %s %s %s %s %s" % ( 
#                                                                                                              ProcessGTExData.SWEAVE_R_SCRIPT, 
#                                                                                                              ProcessGTExData.WORKING_DIR, 
#                                                                                                              ProcessGTExData.ANNOTATION_OUTPUT_FILE, 
#                                                                                                              ProcessGTExData.EXPRESSION_OUTPUT_FILE,
#                                                                                                              ProcessGTExData.EXPRESSION_SAMPLE_OUTPUT_FILE,
#                                                                                                              ProcessGTExData.TX_EXPRESSION_OUTPUT_FOLDER,
#                                                                                                              ProcessGTExData.TX_EXPRESSION_AVG_OUTPUT_FILE
#                                                                                                              )
# 
#         Logger.get_instance().info( "run_statistics : Running command : "+command)
# 
#         self.run_command(command)
