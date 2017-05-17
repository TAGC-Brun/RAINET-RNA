import os
import argparse
import scipy.stats as stats

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer


#===============================================================================
# Started 12-May-2017 
# Diogo Ribeiro
# Based on GroupOddsRatio.py
DESC_COMMENT = "Script to calculate odds ratio of overlap between datasets."
SCRIPT_NAME = "GroupOddsRatio.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read lists of groups from our in-house analysis
# 2) Read lists of groups from external sources
# 3) Read the background (for the statistics)
# 4) Intersect both input lists against background
# 5) Intersect the previous two intersections
# 6) Output odds ratio and pvalues 
#===============================================================================

#===============================================================================
# Processing notes:
# 1) 
#===============================================================================


class GroupOddsRatio(object):

    #=======================================================================
    # Constants
    #=======================================================================
    
    ANNOTATION_FILE_ID_COLUMN = 0
    ANNOTATION_FILE_ANNOTATION_COLUMN = 1
           
    def __init__(self, annotationFile, externalFiles, backgroundList, outputFile):

        self.annotationFile = annotationFile
        self.externalFiles = externalFiles
        self.backgroundList = backgroundList
        self.outputFile = outputFile
                        
        # make output folder
        baseFolder = os.path.basename( self.outputFile)
        if not os.path.exists( baseFolder):
            os.system( baseFolder)


    # #
    # Read list of annotation per RNA.
    def read_annotation_file( self):

        #=======================================================================
        # Example file
        #
        # ENSG00000256751 Predicted
        # ENSG00000256750 Predicted
        # ENSG00000261773 Interacting
        # ENSG00000237402 Interacting
        #=======================================================================
        # The same item can have several annotations

        #=======================================================================
        # initialising
        #=======================================================================

        itemAnnotation = {} # Key -> item ID, value -> set of annotations
        groupItems = {} # Key -> annotation, value -> set of items

        lineCounter = 0

        #=======================================================================
        # read file
        #=======================================================================

        with open( self.annotationFile, "r") as inFile:
            
            for line in inFile:
                line = line.strip()
                lineCounter+=1
                
                spl = line.split( "\t")
                
                ID = spl[ GroupOddsRatio.ANNOTATION_FILE_ID_COLUMN]
                
                # select column to use as annotation                
                annotationItem = spl[ GroupOddsRatio.ANNOTATION_FILE_ANNOTATION_COLUMN]
                    
                # storing tx as key
                if ID not in itemAnnotation:
                    itemAnnotation[ ID] = set()
                itemAnnotation[ ID].add( annotationItem)
                
                # storing annotation as key
                if annotationItem not in groupItems:
                    groupItems[ annotationItem] = set()
                groupItems[ annotationItem].add( ID)
 
            Logger.get_instance().info("GroupOddsRatio.read_annotation_file: number of entries read: %s" % lineCounter )
            Logger.get_instance().info("GroupOddsRatio.read_annotation_file: number of item with annotation: %s" % len( itemAnnotation) )
            Logger.get_instance().info("GroupOddsRatio.read_annotation_file: number of annotations: %s" % len( groupItems) )
         
        for group in sorted(groupItems):
            Logger.get_instance().info("GroupOddsRatio.read_annotation_file %s: %s" % (group, len( groupItems[ group]) ) )
    
 
        self.itemAnnotation = itemAnnotation
        self.groupItems = groupItems


    # #
    # Read each of the provided external files
    def read_external_files(self):
        
        
        # Dictionary which will contain lists of items
        externalLists = {} # key -> name of list (from file name), value -> list of items
        
        # Read file containing list of files
    
        with open( self.externalFiles, "r") as inFile:
            for line in inFile:
                line = line.strip()
                
                # parse name of list
                fileName = line.split("/")[-1]
                fileName = fileName.split( ".")[0]
                
                if fileName not in externalLists:
                    externalLists[ fileName] = set()               
                
                # read the actual file with list of items
                with open( line, "r") as listFile:
                    for ID in listFile:
                        ID = ID.strip()

                        externalLists[ fileName].add( ID)


        Logger.get_instance().info("GroupOddsRatio.read_external_files: number of lists read: %s" % len( externalLists) )
        
        for lst in sorted( externalLists):
            Logger.get_instance().info("GroupOddsRatio.read_external_files %s: %s" % (lst, len( externalLists[ lst]) ) )


        self.externalLists = externalLists


    # #
    # Read background list of items, filters out other lists.
    def read_background_list(self):
        
        #=======================================================================
        # Read file with list of items 
        #=======================================================================

        backgroundItems = set() 

        with open( self.backgroundList, "r") as inFile:
            for line in inFile:
                ID = line.strip()
                
                backgroundItems.add( ID)

        Logger.get_instance().info("GroupOddsRatio.read_background_list: number of items in background: %s" % len( backgroundItems) )

        #=======================================================================
        # Filter for overlap with background items
        #=======================================================================
        
        ### annotation items
                
        annotationItems = self.itemAnnotation.keys()        
        filteredAnnotationItems = [ tx for tx in annotationItems if tx in backgroundItems]
        
        #Note: the overlap should be 100%, but to keep the code stable and flexible it is good to confirm it
        Logger.get_instance().info("GroupOddsRatio.read_background_list: annotation items overlap with background: %s out of %s" % ( len( filteredAnnotationItems), len( annotationItems) ) )
        
        ### External item lists

        filteredExternalLists = {} # key -> name of list (from file name), value -> list of items overlapping with background 

        for lst in self.externalLists:
            
            filteredTxs = [ tx for tx in self.externalLists[ lst] if tx in backgroundItems]
            
            filteredExternalLists[ lst] = filteredTxs
            
            Logger.get_instance().info("GroupOddsRatio.read_background_list: %s items overlap with background: %s out of %s" % ( lst, len( filteredExternalLists[ lst]), len( self.externalLists[ lst]) ) )
            

        self.backgroundItems = backgroundItems
        self.filteredExternalLists = filteredExternalLists
        self.filteredAnnotationItems = filteredAnnotationItems


    # #
    # Calculate odds ratio between groups of lncRNAs, write output file
    def calculate_odds_ratio(self):
        
        #=======================================================================
        # Initialise output files
        #=======================================================================
        
        # output file, format for creating one-side odds ratio table
        outFile = open( self.outputFile, "w")
        outFile.write( "ExternalList\tSize\tItemGroup\tSize\tOverlap\tOddsRatio\tPvalue\n")
        
        # output file, format for creating two-side odds ratio table (yes, no), meaning external list in group, or not.
        outFileTwo = open( self.outputFile + "_two_sided.tsv", "w")
        outFileTwo.write( "ExternalList\tItemGroup\tInGroup\tOverlap\tOddsRatio\tPvalue\n")

        #=======================================================================
        # Produce statistics for each external list
        #=======================================================================

        self.dataStore = {} # stores data for unittesting

        for lst in sorted( self.filteredExternalLists):
            externalList = self.filteredExternalLists[ lst]
            
            # Get overlap between external list and annotated items
            overlapAnnotationExternal = [tx for tx in externalList if tx in self.filteredAnnotationItems]

            #=======================================================================
            # For each annotation group
            #=======================================================================
            for group in sorted( self.groupItems):
                
                # get all the filtered annotated items that belong to current group of annotation
                filteredgroupItems = [tx for tx in self.groupItems[ group] if tx in self.filteredAnnotationItems]

                ## fisher exact test
                # Contigency table:
                #         group    non-group
                # externalList    groupOverlap    len(externalList)-groupOverlap
                # non-externalList    len(group)-groupOverlap    len(backgroundItems) - union(3 other cells)

                # get overlap between annotation and external lists, in the current group
                groupOverlap = {tx for tx in filteredgroupItems if tx in overlapAnnotationExternal}
                
                # get other values for the fisher exact test
                externalNonGroup = set(externalList) - groupOverlap
                nonExternalGroup = set( filteredgroupItems) - groupOverlap
                nonExternalNonGroup = set( self.backgroundItems) - (externalNonGroup.union( nonExternalGroup).union( groupOverlap) )

#                 print lst, group
#                 print len( groupOverlap)
#                 print len( externalNonGroup)
#                 print len( nonExternalGroup)
#                 print len( nonExternalNonGroup)
                
                # assert that sum of columns/rows match
                assert len( externalNonGroup) + len( groupOverlap) == len( externalList)
                assert len( externalNonGroup) + len( nonExternalNonGroup) == len( self.backgroundItems) - ( len( groupOverlap) + len( nonExternalGroup))
                assert len( groupOverlap) + len( nonExternalGroup) == len( filteredgroupItems)
                assert len( groupOverlap) + len( externalNonGroup) + len( nonExternalGroup) + len( nonExternalNonGroup) == len( self.backgroundItems)
    
                ## run the fisher exact test
                matrix = [[len( groupOverlap), len( externalNonGroup)], [len( nonExternalGroup), len( nonExternalNonGroup)]]

                oddsRatio, pvalue = self.fisher_exact_test(matrix)

                self.dataStore[ lst + "|" + group] = matrix

#                 ### Output file in melted format                
#                 # write odds ratio
#                 outFile.write("%s\t%s\t%s\todds_ratio\t%.2f\n" % ( lst, group, len( groupOverlap), oddsRatio))
#                 # write pvalue
#                 outFile.write("%s\t%s\t%s\tpvalue\t%.1e\n" % ( lst, group, len( groupOverlap), pvalue))

                outFile.write("%s\t%s\t%s\t%s\t%s\t%.2f\t%.1e\n" % ( lst, len( externalList), group, len( filteredgroupItems), len( groupOverlap), oddsRatio, pvalue))

                ## output file for two-sided odds ratio table

                # do the opposite test
                matrixInverse = [[len( nonExternalGroup), len( nonExternalNonGroup)], [len( groupOverlap), len( externalNonGroup)] ]
                oddsRatioInverse, pvalueInverse = self.fisher_exact_test(matrixInverse)
                                                
                # write line for odds ratio one way, external list in item group
                outFileTwo.write("%s\t%s\tYes\t%s\t%.2f\t%.1e\n" % ( lst, group, len( groupOverlap), oddsRatio, pvalue))
                
                # write line for odds ratio in the other sense, external list not in item group
                outFileTwo.write("%s\t%s\tNo\t%s\t%.2f\t%.1e\n" % ( lst, group, len( externalNonGroup), oddsRatioInverse, pvalueInverse))
                

        outFile.close() 
        outFileTwo.close()


    # #
    # Perform fishers exact test using scipy
    @staticmethod
    def fisher_exact_test( matrix, alt = "greater"):

        # equivalent of doing in R:
        # fisher.test( matrix(c(8, 2, 1, 5), nrow = 2), alternative = "greater")
        # However: The calculated odds ratio is different from the one R uses.
        # This scipy implementation returns the (more common) unconditional Maximum Likelihood Estimate
        # while R uses the conditional Maximum Likelihood Estimate
    
        oddsRatio, pvalue = stats.fisher_exact( matrix, alternative = alt)
            
        return oddsRatio, pvalue


    # #
    # Run functions in correct order
    def run(self):
        
        Timer.get_instance().step( "Reading annotation file..") 
        run.read_annotation_file( )      

        Timer.get_instance().step( "Reading external files..") 
        run.read_external_files()
        
        Timer.get_instance().step( "Reading background file..") 
        run.read_background_list()

        Timer.get_instance().step( "Creating output file..") 
        run.calculate_odds_ratio()



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
        parser.add_argument('annotationFile', metavar='annotationFile', type=str,
                             help='TSV file with annotation per item. No header. Can have several annotations for same item, one per line. E.g. transcriptID\tannotation.')
        parser.add_argument('externalFiles', metavar='externalFiles', type=str,
                             help='File with list of input files, each with list of items to be used to cross our data. Full paths, no header.')
        parser.add_argument('backgroundList', metavar='backgroundList', type=str,
                             help='File with list of items to use has background. Both annotationFile and externalFiles items will be filtered by this list. No header.')
        parser.add_argument('outputFile', metavar='outputFile', type=str, help='Folder and file where to write output files. Provide full paths.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        # init
        run = GroupOddsRatio( args.annotationFile, args.externalFiles, args.backgroundList, args.outputFile)

        # run functions in order
        run.run()
       
        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

