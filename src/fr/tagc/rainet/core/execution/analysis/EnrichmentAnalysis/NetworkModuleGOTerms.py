import sys
import os
import argparse
import pandas as pd
import scipy.stats as stats

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.NetworkModuleAnnotation import NetworkModuleAnnotation


#===============================================================================
# Started 08-Feb-2017 
# Diogo Ribeiro
# Based on LncRNAGroupAnalysis.py
DESC_COMMENT = "Script to analyse GO terms of enriched network modules."
SCRIPT_NAME = "NetworkModuleGOTerms.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read filtered enrichment results
# 2) Read RAINET DB network modules / GO term correspondence table
# 3) Produce statistics on GO terms 
#===============================================================================

#===============================================================================
# Processing notes:
# 1) 
#===============================================================================


class NetworkModuleGOTerms(object):

    #=======================================================================
    # Constants
    #=======================================================================
    #
           
    def __init__(self, enrichmentFile, rainetDB, outputFolder):

        self.enrichmentFile = enrichmentFile
        self.rainetDB = rainetDB
        self.outputFolder = outputFolder
        
        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()
                
        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)

    # #
    # Get correspondence between transcript IDs and gene IDs using rainetDB
    def read_rainet_db(self):

        moduleGOTerms = {} # key -> network module ID, val -> set of associated GO terms
        goTermFrequency = {} # key -> GO term, val -> set of network module IDs
        termDescriptions = {} # key -> GO term, val -> description of GO term
 
        #===============================================================================
        # Query NetworkMolduleAnnotation table, which contains GO terms associated to each NetworkModule
        #===============================================================================
        query = self.sql_session.query( NetworkModuleAnnotation.networkModule_id, NetworkModuleAnnotation.annotationID, NetworkModuleAnnotation.annotationTerm ).all()
        # Note: a networkModule can have several GO terms, a GO term may occur in many networkModules      
                
        for moduleID, goTerm, termDescription in query:
            
#             # skip entries that are empty
#             if goTerm == "":
#                 continue
            
            moduleID = str( moduleID)
            goTerm = str( goTerm)
            termDescription = str( termDescription)
            
            if moduleID not in moduleGOTerms:
                moduleGOTerms[ moduleID] = set()
            moduleGOTerms[ moduleID].add( goTerm)
            
            if goTerm not in goTermFrequency:
                goTermFrequency[ goTerm] = set()
            goTermFrequency[ goTerm].add( moduleID)

            # associate GO term ID to a description            
            termDescriptions[ goTerm] = termDescription
        
        assert( len( termDescriptions) == len( goTermFrequency))

        Logger.get_instance().info( "NetworkModuleGOTerms.read_rainet_db : Number of Module-GOterm associations in DB: %s" % ( len( query)) )
        Logger.get_instance().info( "NetworkModuleGOTerms.read_rainet_db : Number Modules with annotation %s" % ( len( moduleGOTerms)) )
        Logger.get_instance().info( "NetworkModuleGOTerms.read_rainet_db : Number GO terms associated to module %s" % ( len( goTermFrequency)) )
     
        self.moduleGOTerms = moduleGOTerms
        self.goTermFrequency = goTermFrequency
        self.termDescriptions = termDescriptions


    # #
    # Read enrichment file and produce output file
    def read_enrichment_file(self):
        
        #=======================================================================
        # Enrichment file
        #=======================================================================
        # example format: 
        #
        # transcriptID    annotID number_observed_interactions    number_possible_interactions    percent_interacting_proteins    warning pval    corrected_pval  sign_corrected
        # ENST00000354541 235     5       25      4.6%    0       7.7e-04 4.4e-02 1
        # ENST00000354541 376     10      73      4.6%    0       4.6e-04 3.3e-02 1

        # the important information is the annotID, which indicates which NetworkModule is enriched.
        # Note that the same network modules will show up for several transcripts, and this gives us a frequency of GO term presence
        
        
        #=======================================================================
        # output file
        #=======================================================================
        # example format: 
        # GO_term\tterm_description\tenrichment_occurrences\tmodules_with_go_term\toccurrence_per_module\n
        
        outFile = open( self.outputFolder + "/module_GO_term_occurrence.tsv", "w")
        
        outFile.write( "GO_term\tterm_description\tenrichment_occurrences\tmodules_with_go_term\toccurrence_per_module\n")
        
        # one line for each GO term with data
        
        termOccurence = {} # Key -> GO term, val -> how many times this GO term appears on enrichment file
        
        nLines = 0
        with open( self.enrichmentFile, "r") as inFile:
            
            header = inFile.readline()
            
            for line in inFile:
                spl = line.strip().split( "\t")
                
                nLines+=1
                
                moduleID = spl[ 1]
                
                if moduleID in self.moduleGOTerms:
                    for term in self.moduleGOTerms[ moduleID]:
                        #print term, self.termDescriptions[ term]
                        if term not in termOccurence:
                            termOccurence[ term] = 0
                            
                        termOccurence[ term] += 1

                else:
                    # modules that have no GO term annotations
                    # does not happen if allowing processing of empty GO terms
                    pass
        

        Logger.get_instance().info( "NetworkModuleGOTerms.read_enrichment_file : read %s enrichments." % ( nLines) )
        Logger.get_instance().info( "NetworkModuleGOTerms.read_enrichment_file : %s GO terms with data." % ( len( termOccurence) ) )
        
        # sort occurrence results by ratio of occurrence in enrichment file versus occurrence in the population of modules        
        sortedByRatio = sorted( termOccurence.items(), key=lambda x: float( x[1]) / len( self.goTermFrequency[ x[0]]), reverse = True)
        
        for term, value in sortedByRatio:
            ratio = float(value)/len( self.goTermFrequency[ term])
            outFile.write( "%s\t%s\t%i\t%i\t%.2f\n" % (term, self.termDescriptions[ term], value, len( self.goTermFrequency[ term]), ratio) )

        outFile.close()

        # awk '$2=="840"' enrichment_results_filtered.tsv | wc -l  --> 66


#     # #
#     # Calculate odds ratio between groups of lncRNAs, write output file
#     def calculate_odds_ratio(self):
#         
#         #=======================================================================
#         # Output file
#         #=======================================================================
#         # format: melted file 
#         # ExternalList\tTranscriptGroup\tMetric\tValue
#         outFile = open( self.outputFile, "w")
# #        outFile.write( "ExternalList\tTranscriptGroup\tOverlap\tMetric\tValue\n")
#         outFile.write( "ExternalList\tTranscriptGroup\tOverlap\tOddsRatio\tPvalue\n")
# 
#         #=======================================================================
#         # Produce statistics for each external list
#         #=======================================================================
# 
#         self.dataStore = {} # stores data for unittesting
# 
#         for lst in sorted( self.filteredExternalLists):
#             externalList = self.filteredExternalLists[ lst]
#             
#             # Get overlap between external list and annotated transcripts
#             overlapAnnotationExternal = [tx for tx in externalList if tx in self.filteredAnnotationTranscripts]
# 
#             #=======================================================================
#             # For each annotation group
#             #=======================================================================
#             for group in sorted( self.groupTranscripts):
#                 
#                 # get all the filtered annotated transcripts that belong to current group of annotation
#                 filteredGroupTranscripts = [tx for tx in self.groupTranscripts[ group] if tx in self.filteredAnnotationTranscripts]
# 
#                 ## fisher exact test
#                 # Our question: does the current annotation genes (e.g. interacting, enriched) significantly overlap with the genes of the current external list (Liu2016, Mukherjee2016)?
#                 # Contigency table:
#                 #         group    non-group
#                 # externalList    groupOverlap    len(externalList)-groupOverlap
#                 # non-externalList    len(group)-groupOverlap    len(backgroundTranscripts) - union(3 other cells)
# 
#                 # get overlap between annotation and external lists, in the current group
#                 groupOverlap = {tx for tx in filteredGroupTranscripts if tx in overlapAnnotationExternal}
#                 
#                 # get other values for the fisher exact test
#                 externalNonGroup = set(externalList) - groupOverlap
#                 nonExternalGroup = set( filteredGroupTranscripts) - groupOverlap
#                 nonExternalNonGroup = set( self.backgroundTranscripts) - (externalNonGroup.union( nonExternalGroup).union( groupOverlap) )
# 
# #                 print lst, group
# #                 print len( groupOverlap)
# #                 print len( externalNonGroup)
# #                 print len( nonExternalGroup)
# #                 print len( nonExternalNonGroup)
#                 
#                 # assert that sum of columns/rows match
#                 assert len( externalNonGroup) + len( groupOverlap) == len( externalList)
#                 assert len( externalNonGroup) + len( nonExternalNonGroup) == len( self.backgroundTranscripts) - ( len( groupOverlap) + len( nonExternalGroup))
#                 assert len( groupOverlap) + len( nonExternalGroup) == len( filteredGroupTranscripts)
#                 assert len( groupOverlap) + len( externalNonGroup) + len( nonExternalGroup) + len( nonExternalNonGroup) == len( self.backgroundTranscripts)
#     
#                 # run the fisher exact test
#                 matrix = [[len( groupOverlap), len( externalNonGroup)], [len( nonExternalGroup), len( nonExternalNonGroup)]]
# 
#                 oddsRatio, pvalue = self.fisher_exact_test(matrix)
# 
#                 self.dataStore[ lst + "|" + group] = matrix
# 
# #                 ### Output file in melted format                
# #                 # write odds ratio
# #                 outFile.write("%s\t%s\t%s\todds_ratio\t%.2f\n" % ( lst, group, len( groupOverlap), oddsRatio))
# #                 # write pvalue
# #                 outFile.write("%s\t%s\t%s\tpvalue\t%.1e\n" % ( lst, group, len( groupOverlap), pvalue))
# 
#                 outFile.write("%s\t%s\t%s\t%.2f\t%.1e\n" % ( lst, group, len( groupOverlap), oddsRatio, pvalue))
#                 
# 
#         outFile.close() 
# 
# 
#     # #
#     # Perform fishers exact test using scipy
#     @staticmethod
#     def fisher_exact_test( matrix):
# 
#         # equivalent of doing in R:
#         # fisher.test( matrix(c(8, 2, 1, 5), nrow = 2), alternative = "greater")
#         # However: The calculated odds ratio is different from the one R uses.
#         # This scipy implementation returns the (more common) unconditional Maximum Likelihood Estimate
#         # while R uses the conditional Maximum Likelihood Estimate
#     
#         oddsRatio, pvalue = stats.fisher_exact( matrix, alternative = "greater")
#             
#         return oddsRatio, pvalue


    # #
    # Run functions in correct order
    def run(self):
        
        Timer.get_instance().step( "Reading RAINET DB..") 
        self.read_rainet_db( )      

        Timer.get_instance().step( "Reading enrichment file..") 
        self.read_enrichment_file()


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
        parser.add_argument('enrichmentFile', metavar='enrichmentFile', type=str,
                             help='Rainet database to be used for ID mapping.')       
        parser.add_argument('rainetDB', metavar='rainetDB', type=str,
                             help='Rainet database to be used for ID mapping.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        # init
        run = NetworkModuleGOTerms( args.enrichmentFile, args.rainetDB, args.outputFolder)

        # run functions in order
        run.run()
       
        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

