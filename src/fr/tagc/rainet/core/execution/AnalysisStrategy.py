
from sqlalchemy import or_,and_,distinct

from fr.tagc.rainet.core.execution.ExecutionStrategy import ExecutionStrategy
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.util import Constants as OptionConstantsConstants
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


# #
# This class define the Strategy enabling user to perform DB query interactively
class AnalysisStrategy(ExecutionStrategy):

    # #
    # The Strategy execution method
    def execute(self):

        # Getting input arguments        
        self.DBPath = OptionManager.get_instance().get_option(OptionConstants.OPTION_DB_NAME)
        self.species = OptionManager.get_instance().get_option(OptionConstants.OPTION_SPECIES)
        self.minimumInteractionScore =  OptionManager.get_instance().get_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE)
        self.transcriptBiotype = OptionManager.get_instance().get_option(OptionConstants.OPTION_TRANSCRIPT_BIOTYPE)
        self.lncRNABiotypes = OptionManager.get_instance().get_option(OptionConstants.OPTION_LNCRNA_BIOTYPES)
        self.gencode = OptionManager.get_instance().get_option(OptionConstants.OPTION_GENCODE)
        
        # Check if minimumInteractionScore is float or OFF
        if self.minimumInteractionScore != OptionConstantsConstants.DEFAULT_INTERACTION_SCORE:
            try:
                float(self.minimumInteractionScore)
            except TypeError:
                raise RainetException( "AnalysisStrategy.execute: Provided minimum interaction score is not a float.")

        # Check if transcript biotype is correct
        if self.transcriptBiotype not in OptionConstantsConstants.TRANSCRIPT_BIOTYPES:
            raise RainetException( "AnalysisStrategy.execute: Provided transcript biotype is not allowed: "+str(self.transcriptBiotype))

        # Process and check provided list of lncRNA subtypes
        if self.lncRNABiotypes != OptionConstantsConstants.DEFAULT_LNCRNA_BIOTYPES:
            bioTypeStr = ""
            for subtype in self.lncRNABiotypes.split(","):
                if subtype not in OptionConstantsConstants.LNCRNA_BIOTYPES:
                    raise RainetException( "AnalysisStrategy.execute: Provided lncRNA biotype is not allowed: "+str(subtype))              
                bioTypeStr+='"'+subtype+'",'
            self.lncRNABiotypes = bioTypeStr[:-1] # taking out last comma

        # Check if gencode argument is correct
        if self.gencode != OptionConstantsConstants.DEFAULT_GENCODE:
            try:
                self.gencode = int(self.gencode)
            except TypeError:
                raise RainetException( "AnalysisStrategy.execute: Provided GencodeBasicOnly argument is not numeric.")
            if self.gencode != 1:
                raise RainetException( "AnalysisStrategy.execute: Provided GencodeBasicOnly argument must be either OFF or 1.")

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.DBPath)
        self.sql_session = SQLManager.get_instance().get_session()
        self.db_engine = SQLManager.get_instance().get_engine()
                        
        self.analysis()
        
    # #
    # Run analysis-related functions in order
    def analysis(self):

        Logger.get_instance().info( "AnalysisStrategy.analysis: Starting..." )
        
        print ("ARGUMENTS:\n"+"\n".join([self.DBPath,self.species,self.minimumInteractionScore,self.transcriptBiotype,self.lncRNABiotypes]) )

        #===================================================================
        # Filter datasets
        #===================================================================

        self.filterRNA()
        
        self.filterProtein()
        
        self.filterPRI()

        #===================================================================
        # Perform analysis
        #===================================================================

        self.enrichementAnalysis()


    # #
    # Filter RNA models
    #
    # Stores final list of RNAs on DataManager 
    def filterRNA(self):

        # Get all RNA items
        DataManager.get_instance().perform_query("allRNAs", "query(RNA).all()")
        allRNAs = DataManager.get_instance().get_data("allRNAs")
        
        # Filter transcripts based on biotype
        if self.transcriptBiotype != OptionConstantsConstants.BIOTYPE_LNCRNA:        
            query = "query("+self.transcriptBiotype+").all()"
        elif self.lncRNABiotypes == OptionConstantsConstants.DEFAULT_LNCRNA_BIOTYPES:
            query = "query("+self.transcriptBiotype+").all()"
        else:
        # Filter LncRNAs based on wanted subtypes (if lncRNA biotype chosen)
            query = "query("+self.transcriptBiotype+").filter("+self.transcriptBiotype+".transcriptBiotype.in_(["+self.lncRNABiotypes+"])).all() "

        DataManager.get_instance().perform_query("filterRNA1", query) 
        filterRNA1 = {str(item.transcriptID) for item in DataManager.get_instance().get_data("filterRNA1")}

        # Filter transcripts based on gencode_basic presence
        if self.gencode == 1:
            query = "query(RNA).filter(RNA.transcriptGencodeBasic == 1).all()"
        else:
            query = "query(RNA).all()"
        
        DataManager.get_instance().perform_query("filterRNA2", query)  
        filterRNA2 = {str(item.transcriptID) for item in DataManager.get_instance().get_data("filterRNA2")}
                       
        # Select one transcript isoform per gene
        # TODO

        # Get intersection of the various filterings
        selectedRNAs = []
        for rna in allRNAs:
            if rna.transcriptID in filterRNA1 and rna.transcriptID in filterRNA2:
                selectedRNAs.append(rna)

        DataManager.get_instance().store_data(DataConstants.RNA_FILTER_KW, selectedRNAs)


    # #
    # Filter protein models    
    #
    # Stores final list of proteins on DataManager 
    def filterProtein(self):

        query = "query( Protein ).all()"
        DataManager.get_instance().perform_query(DataConstants.PROT_FILTER_KW, query) 

        # Will filter out:
        # - Protein isoforms (?)
        # - Proteins based on their mRNA expression


    # #
    # Filter protein-RNA interactions
    def filterPRI(self):

        # Retrieve selected RNAs and proteins
         
        selectedRNAs = {str(item.transcriptID) for item in DataManager.get_instance().get_data(DataConstants.RNA_FILTER_KW) }
        selectedProteins = {str(item.uniprotAC) for item in DataManager.get_instance().get_data(DataConstants.PROT_FILTER_KW) } 
         
        # Filter interactions based on minimumInteractionScore
         
        if self.minimumInteractionScore != OptionConstantsConstants.DEFAULT_INTERACTION_SCORE:
            query = "query( ProteinRNAInteractionCatRAPID ).filter(ProteinRNAInteractionCatRAPID.interactionScore >= "+str(self.minimumInteractionScore)+").all()"    
        else:
            query = "query( ProteinRNAInteractionCatRAPID ).all()"
 
        DataManager.get_instance().perform_query(DataConstants.PRI_FILTER_KW, query)
 
        interactions = DataManager.get_instance().get_data(DataConstants.PRI_FILTER_KW)
 
        # Filter for interactions between selected RNAs and proteins
        selectedInteractions = []
        for inter in interactions:
            if inter.transcriptID in selectedRNAs and inter.proteinID in selectedProteins:
                selectedInteractions.append(inter)
 
        DataManager.get_instance().store_data(DataConstants.PRI_FILTER_KW, selectedInteractions)

        
    def enrichementAnalysis(self):

        pass

# 
#         # For each interacting RNA, does it target significantly more proteins in same KEGG pathway
#         
#         # 1) for each KEGG pathway, count how many proteins in it
#         keggFrequencies = {}
# 
#         keggPathways = self.sql_session.query( KEGGPathway ).all()
#         for pathway in keggPathways:
#             print (pathway)

#         selectedInteractions = DataManager.get_instance().get_data(self.PRIKeyword)
#         # Get protein IDs and transcript IDs present in filtered interactions
#         for inter in selectedInteractions:
#             print (inter.transcriptID,inter.proteinID)
#              
#             prot = self.sql_session.query( Protein ).filter( Protein.uniprotAC == inter.proteinID ).all()
#  
#             rna = self.sql_session.query( RNA ).filter( RNA.transcriptID == inter.transcriptID ).all()


