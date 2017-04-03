import sys
import os
import argparse

import numpy as np
# import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

# from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
# from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
# from fr.tagc.rainet.core.util.data.DataManager import DataManager

# from fr.tagc.rainet.core.data.Protein import Protein

#===============================================================================
# Started 25-June-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to parse output from EnrichmentAnalysisStrategy and add extra filters as well as transform into matrix format."
SCRIPT_NAME = "FilterEnrichmentResults.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read enrichment_per_rna file, get list of RNAs to filter out based on input criteria.
# 2) Rewrite enrichment_results file based on these filters, output also enrichment matrix file.
# 3) Options to apply other types of filtering
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Only writing to file annotations with at least one enrichment
# 2) If random has zero significant, ratio is close to infinite
# 3) When using topEnrichmentsPerComplex option, get all enrichments with same number of observed interactions in case of draw
# 4) Output matrix file is written after the enrichment_per_rna random filtering and minimumProteinInteraction filtering, but before other optional filters
#===============================================================================


class FilterEnrichmentResults(object):

    #===============================================================================
    # Class constants
    #===============================================================================
    
    REPORT_LIST_RNA_SIGN_ENRICH = "/list_RNA_above_random.txt"
    REPORT_FILTERED_RNA_ANNOT_RESULTS = "/enrichment_results_filtered.tsv"
    REPORT_RNA_ANNOT_RESULTS_MATRIX = "/enrichment_results_filtered_matrix.tsv"
    REPORT_MATRIX_ROW_ANNOTATION = "/matrix_row_annotation.tsv"
    REPORT_MATRIX_COL_ANNOTATION = "/matrix_col_annotation.tsv"
    REPORT_SPECIFICITY_RANK = "/enrichment_specificity_rank.tsv"
    REPORT_ENRICHMENT_SUMMARY = "/enrichment_summary.tsv"
    
    SEVERAL_ANNOTATION_TAG = "Overlapping_annotations"
    
    COLORS_SET3 = ["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999"]
    COLORS_SET2 = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"]
        
    def __init__(self, enrichmentPerRNAFile, enrichmentResultsFile, outputFolder, matrixValueColumn, filterWarningColumn, \
                      filterWarningValue, minimumRatio, rowAnnotationFile, colAnnotationFile, \
                      maskMultiple, noAnnotationTag, noAnnotationFilter, annotSpecificityFilter, \
                      transcriptSpecificityFilter, minimumProteinInteraction, topEnrichmentsPerComplex):

        self.enrichmentPerRNAFile = enrichmentPerRNAFile
        self.enrichmentResultsFile = enrichmentResultsFile
        self.outputFolder = outputFolder
        self.matrixValueColumn = matrixValueColumn
        self.filterWarningColumn = filterWarningColumn
        self.filterWarningValue = filterWarningValue
        self.minimumRatio = minimumRatio
        self.rowAnnotationFile = rowAnnotationFile
        self.colAnnotationFile = colAnnotationFile
        self.maskMultiple = maskMultiple
        self.noAnnotationTag = noAnnotationTag
        self.noAnnotationFilter = noAnnotationFilter
        self.annotSpecificityFilter = annotSpecificityFilter
        self.transcriptSpecificityFilter = transcriptSpecificityFilter
        self.minimumProteinInteraction = minimumProteinInteraction
        self.topEnrichmentsPerComplex = topEnrichmentsPerComplex

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)
            
        # other non input variables
        self.colAnnotation = {}
        self.rowAnnotation = {}
    
    # #
    # Read RNA enrichment file, get list of RNAs with significantly more enrichments compared to control.
    def read_enrichment_per_rna_file( self):
           
        # Example format:
        # transcriptID    n_sign_tests_no_warning avg_n_sign_random_no_warning    n_times_above_random    empiricalPvalue significant
        # ENST00000230113 57      50.78   815     1.9e-01 0
    
        poolRNA = set()
    
        listRNASignificantEnrich = set()
    
        # before filtering
        observedValues = []
        randomValues = []
        diffValues = []
    
        # after filtering
        observedValuesAfter = []
        randomValuesAfter = []
    
        countAboveRandom = 0
        countRatioPassed = 0
        
        with open( self.enrichmentPerRNAFile, "r") as inFile:
            
            inFile.readline() # skip header
    
            for line in inFile:
                line = line.strip()
    
                # unpack line            
                txID, observedSign, randomSign, _, _, signFlag = line.split( "\t")
                
                poolRNA.add( txID)
    
                # booleans for filterings, start as positive
                signBoo = 1
                ratioBoo = 1           
    
                if signFlag == "0":
                    signBoo = 0
                else:
                    countAboveRandom += 1
    
                if self.minimumRatio != "OFF":
                    # filter also by ratio between real and random
                    randomSign = float( randomSign)
                    
                    if randomSign != 0: # to avoid division by Zero
                        ratio = float( observedSign) / randomSign
                    else:
                        ratio = float( observedSign) / 0.000000000000000000000000000000000000000001
    
                    if ratio < self.minimumRatio:
                        ratioBoo = 0
                    else:
                        countRatioPassed += 1                    
    
                # if both filters passed
                if signBoo and ratioBoo:
                    listRNASignificantEnrich.add( txID)
                    observedValuesAfter.append( float(observedSign))
                    randomValuesAfter.append( float(randomSign))
                    
    
                observedValues.append( float(observedSign))
                randomValues.append( float(randomSign))
    
                # difference between observed and random number of tests with enrichment
                diffValue = float( observedSign) - float( randomSign)
                diffValues.append( diffValue)
    
        # Print basic stats
    
        countRealEnrichments = sum(observedValues)
        countRandomEnrichments = sum(randomValues)
        
        Logger.get_instance().info( "read_enrichment_per_rna_file : Total number of RNAs: %s " % len( poolRNA) )
    
        Logger.get_instance().info( "read_enrichment_per_rna_file : Number of RNAs with enrichment significantly above random: %s " % countAboveRandom )
        if self.minimumRatio != "OFF":
            Logger.get_instance().info( "read_enrichment_per_rna_file : Number of RNAs passing real/random ratio: %s " % countRatioPassed )
        Logger.get_instance().info( "read_enrichment_per_rna_file : Number of RNAs passing all filters: %s " % len( listRNASignificantEnrich) )
    
        Logger.get_instance().info( "read_enrichment_per_rna_file : Observed enrichments: %i" % ( countRealEnrichments ) )
        Logger.get_instance().info( "read_enrichment_per_rna_file : Random enrichments: %.1f" % ( countRandomEnrichments) )
        
        Logger.get_instance().info( "read_enrichment_per_rna_file : Significant enrichment tests:" )
        Logger.get_instance().info( "read_enrichment_per_rna_file : (before filtering) Observed mean: %.1f Random mean: %.1f Mean difference: %.1f" % ( np.mean( observedValues), np.mean( randomValues), np.mean( diffValues)) )
        Logger.get_instance().info( "read_enrichment_per_rna_file : (after filtering) Observed mean: %.1f Random mean: %.1f" % ( np.mean( observedValuesAfter), np.mean( randomValuesAfter)) )
        
        return listRNASignificantEnrich, countRealEnrichments, countRandomEnrichments

        
    # #
    # Read RNA-annotation enrichment file (results), filter out results for RNAs that are not significantly enriched over random control, produce output files.
    # Optional, further filtering of enrichments with warning flag on, and minimum of protein interactions
    def read_enrichment_results_file(self, list_rna_significant_enrich, count_real_enrichments, count_random_enrichments):
        
        # Example format:
        # transcriptID    annotID number_observed_interactions    number_possible_interactions    total_interacting_proteins      warning pval    corrected_pval  sign_corrected
        # ENST00000230113 344     12      15      7515    0       5.0e-01 7.3e-01 0   
    
        #===============================================================================
        # Output files
        #===============================================================================
    
        # File with enrichment_results file data in matrix format.
        # only for RNAs that pass enrichment_per_rna filter
        outFile2 = open( self.outputFolder + FilterEnrichmentResults.REPORT_RNA_ANNOT_RESULTS_MATRIX, "w")
    
        #===============================================================================
        # Read file, apply filter, write output
        #===============================================================================
    
        # initialise some variables
        excludedByRNA = 0
        excludedByMinimumInteractions = 0
        setRNAs = set()
        setAnnots = set()
        dictPairs = {} # key -> txID|annotID, val -> value (variable, boolean or pval) # includes annotations without enrichments
        annotValues = {} # key -> annot ID, val -> list of values # includes annotations without enrichments
    
        filteredEnrichmentResults = [] # lists with enrichment results after filter (same as output file)
    
        with open( self.enrichmentResultsFile, "r") as inFile:
            
            header = inFile.readline()
#             outFile1.write( header) # transport header
            filteredEnrichmentResults.append( header)
    
            for line in inFile:
                
                # First check if all the results with this RNA should be filtered out or not
                txID = line.split("\t")[0]
    
                if txID not in list_rna_significant_enrich:
                    excludedByRNA += 1
                    continue
    
                # RNA was not filtered out:    
                line = line.strip()
                spl = line.split( "\t")
                
                annotID = spl[1]
                nObservedInteractions = int( spl[2])
                warningFlag = int( spl[5])
                signFlag = int( spl[8])
                
                value = spl[ self.matrixValueColumn] # can be p-value, corrected p-value or significant boolean
    
                # if filtering is on
                if self.filterWarningColumn:
                    # If value is determined not significant OR is tagged with a warning (for several reasons), fill it with constant value. 
                    if warningFlag or signFlag == 0:
                        value = self.filterWarningValue
                    # apply filter for minimum number of interactions in enrichment
                    elif self.minimumProteinInteraction != -1 and nObservedInteractions < self.minimumProteinInteraction:
                        excludedByMinimumInteractions += 1
                        value = self.filterWarningValue
                    else:
                        filteredEnrichmentResults.append( line)
    
                pair = txID + "|" + annotID
                if pair not in dictPairs:
                    dictPairs[ pair] = value
                else:
                    raise RainetException("read_enrichment_results_file: duplicate key", pair)
    
                if annotID not in annotValues:
                    annotValues[ annotID] = []
                annotValues[ annotID].append( value)
                       
                setRNAs.add( txID)
                setAnnots.add( annotID)
    
    
        # Get set of annotations which had no enrichment
        nonEnrichedAnnotations = set()
        for annotID in annotValues:
            if annotValues[ annotID].count( self.filterWarningValue) == len( annotValues[ annotID]):
                nonEnrichedAnnotations.add( annotID)
    
    
        # the number of retained enrichments after filtering is equal to the length of filteredEnrichmentsResults minus the header
        countFilteredEnrichments = len(filteredEnrichmentResults) - 1
    
        Logger.get_instance().info( "read_enrichment_results_file : Number of lines filtered out because of enrichment_rna filter: %s" % ( excludedByRNA) )
        Logger.get_instance().info( "read_enrichment_results_file : Number of lines filtered out because of minimum observed proteins interacting filter (after previous filter): %s" % ( excludedByMinimumInteractions) )
        Logger.get_instance().info( "read_enrichment_results_file : Total number of enrichments retained after filtering: %s" % ( countFilteredEnrichments) )
    
        assert len( dictPairs) ==  len( setRNAs) * len( setAnnots), "number of pairs should equal number of RNAs times number of annotations"
        
    
        #===============================================================================
        # Writing matrix file and associated files
        #
        #===============================================================================
    
        # initialise
        sortedSetRNAs = sorted( setRNAs)   
        sortedSetAnnots = sorted( setAnnots - nonEnrichedAnnotations)
    
        #===============================================================================
        # Writing row annotation file (if applicable)
        #===============================================================================
    
        if len( self.rowAnnotation) > 0: # if rowAnnotation option is on
            rowWithAnnotation = self._write_matrix_annotation_file( FilterEnrichmentResults.REPORT_MATRIX_ROW_ANNOTATION, sortedSetRNAs, self.rowAnnotation, self.maskMultiple, self.noAnnotationTag, FilterEnrichmentResults.COLORS_SET2, self.noAnnotationFilter)
        else:
            rowWithAnnotation = None
    
        #===============================================================================
        # Writing col annotation file (if applicable)
        #===============================================================================
    
        if len( self.colAnnotation) > 0: # if colAnnotation option is on
            colWithAnnotation = self._write_matrix_annotation_file( FilterEnrichmentResults.REPORT_MATRIX_COL_ANNOTATION, sortedSetAnnots, self.colAnnotation, self.maskMultiple, self.noAnnotationTag, FilterEnrichmentResults.COLORS_SET3, self.noAnnotationFilter)
        else:
            colWithAnnotation = None
    
        #===============================================================================
        # Writing matrix file
        #===============================================================================
    
        # if no_annotation filter is active, filter lists of RNAs and Annotations accordingly to write matrix
        if self.noAnnotationFilter:
            sortedSetRNAs = [ item for item in sortedSetRNAs if item in rowWithAnnotation]
            sortedSetAnnots = [ item for item in sortedSetAnnots if item in colWithAnnotation]
    
        # write header with annot IDs
        outFile2.write( "RNAs")
        for annot in sortedSetAnnots:
            outFile2.write( "\t%s" % annot )
        outFile2.write( "\n")
           
        # write bulk of file, one row per rna, one column per protein
        for rna in sortedSetRNAs:
            text = rna
            for annot in sortedSetAnnots:          
                tag = rna + "|" + annot
                if tag in dictPairs:
                    score = dictPairs[tag]
                else:
                    raise RainetException("read_enrichment_results_file: no data for pair:", pair)
                text+= "\t%s" % score 
            text+= "\n"
            outFile2.write( text)
        
        Logger.get_instance().info( "read_enrichment_results_file : Number of RNAs in matrix file (rows): %i" % ( len( sortedSetRNAs) ) )
        Logger.get_instance().info( "read_enrichment_results_file : Number of Annotations in matrix file (columns): %i" % ( len( sortedSetAnnots) ) )
        
        outFile2.close()
        
        return dictPairs, filteredEnrichmentResults
    
    
    # #
    # Writes file with complex-lncRNA pairs ranked by specificity
    def rank_by_specificity(self, dict_pairs, filter_warning_value):
    
        lncDict = {} # key -> txID, val -> list of annot ID where enriched
        annotDict = {} # key -> annot ID, val -> list of lncRNAs where enriched
    
        # read lncRNA-annot pairs with enrichments, store data to have counts of enrichments for each lncRNA and each annotation    
        for pair in dict_pairs:
            txID, annotID = pair.split("|")
            
            if dict_pairs[ pair] != filter_warning_value:
                        
                if txID not in lncDict:
                    lncDict[ txID] = []
        
                lncDict[ txID].append( annotID)
        
                if annotID not in annotDict:
                    annotDict[ annotID] = []
        
                annotDict[ annotID].append( txID)
    
    
        # read again dictPairs and put results into a rank so that pairs can be ranked by 'global' specificity
        rankDict = {} # key -> rank value of lncRNA-annot specificity, val -> lncRNA-complex pair
        for pair in dict_pairs:
            txID, annotID = pair.split("|")
            if dict_pairs[ pair] != filter_warning_value:
    
                rank = len( lncDict[ txID]) + len( annotDict[ annotID])
    
                if rank not in rankDict:
                    rankDict[ rank] = []
                rankDict[ rank].append( pair)
    
        #===============================================================================
        # Writing specificity rank file
        #===============================================================================
        # e.g. 
        # lncRNA\tcomplex\tlncRNA_specificity\tannot_specificity\n
        # lnc1\ttelomerase\t2\t5\n
    
        # File with same file as enrichment_results file but filtered based on enrichment_per_rna
        outFile = open( self.outputFolder + FilterEnrichmentResults.REPORT_SPECIFICITY_RANK, "w")
    
        outFile.write( "transcriptID\tannotID\ttranscript_enrichments\tannot_enrichments\n")
    
        for rank in sorted( rankDict):
    
            for pair in rankDict[ rank]:
                txID, annotID = pair.split("|")
                
                outFile.write( "%s\t%s\t%s\t%s\n" % ( txID, annotID, len( lncDict[ txID]), len( annotDict[ annotID])) )
    
        outFile.close()
    
        return annotDict, lncDict
    
    
    # #
    # Helper function to write annotation files associated to matrix file.
    # If option chose, filters out items with no annotation
    def _write_matrix_annotation_file(self, output_file, sorted_list, annotation_dict, mask_multiple, no_annotation_tag, color_set, no_annotation_filter):
    
        Logger.get_instance().info( "_write_matrix_annotation_file : writing %s file.." % ( output_file ) )
    
        # set up colors for category plotting
        listOfColors = []
        lookupColors = {}
        colorCount = 0
    
        withAnnotation = set() # stores items with annotation
    
        annotations = [] # contains list of annotations in order
        for item in sorted_list:
            if item in annotation_dict:
                # if there is several annotations and we only want information for proteins with a single domain
                if len( annotation_dict[ item]) > 1 and mask_multiple == 1:
                    annotation = FilterEnrichmentResults.SEVERAL_ANNOTATION_TAG
                # if wanting information for all annotations
                elif len(annotation_dict[ item]) > 0:
                    # get all annotations
                    annotation = ",".join( list( annotation_dict[ item]) )
                else:
                    raise RainetException("_write_matrix_annotation_file: Annotation information is incorrect. ", annotation_dict[ item])
            # if there is no annotation information for current item
            else:
                annotation = no_annotation_tag
    
            # first time annotation is seen, attribute color to it
            if annotation not in lookupColors:
                try:
                    lookupColors[ annotation] = color_set[ colorCount]
                except IndexError:
                    raise RainetException("_write_matrix_annotation_file: Too many different categories for the number of possible colors. ")
                colorCount+= 1
        
            if no_annotation_filter and annotation == no_annotation_tag:
                # if no_annotation filter on and there is no annotation, do not store in variables
                pass
            else:
                annotations.append( annotation)
                listOfColors.append( lookupColors[ annotation])
                withAnnotation.add( item)        
    
            assert len( listOfColors) == len( annotations)
            assert len( set(listOfColors)) == len( set( annotations))
    
        # Write file with annotation of rows (default = RNA) based on same sorting as output matrix (for R script)
        outFile = open( self.outputFolder + output_file, "w")        
        outFile.write( "\t".join( annotations) + "\n")
        outFile.write( "\t".join( listOfColors) + "\n")
        outFile.close()
    
        Logger.get_instance().info( "_write_matrix_annotation_file : Number of items filtered out by having no annotation: %i" % ( len( sorted_list) - len( withAnnotation) ) )
    
        return withAnnotation
    
    # #
    # Reads annotation file, for adding information on columns or rows of enrichment matrix 
    def read_annotation_file(self, annotation_file):
        
        annotationDict = {} # key -> RNA/Protein group identifier, value -> set of annotations
        
        with open( annotation_file, "r") as inFile:
            for line in inFile:
                line = line.strip()
                spl = line.split( "\t")
                ID = spl[0]
                annotText = spl[1]
                
                if ID not in annotationDict:
                    annotationDict[ ID] = set()
                annotationDict[ ID].add( annotText)
    
        Logger.get_instance().info( "read_annotation_file : %s file read. Number of items with annotation: %s " % ( annotation_file , len( annotationDict) ) )
    
        return annotationDict
    
    
    # #
    # Use specificity ranking and list of previously filtered enrichments, filter them further based on annotation specificity levels.
    def filter_by_specificity(self, filtered_enrichment_results, annot_dict, lnc_dict):
         
        notPassingFilters = 0

        filteredEnrichmentResults = []
        
        # write header
        filteredEnrichmentResults.append( filtered_enrichment_results[0])
        
        for enrich in filtered_enrichment_results[1:]:
          
            spl = enrich.split( "\t")
             
            annotID = spl[1]
            transcriptID = spl[0]
             
            # keep only enrichments where number of enriched transcripts is lower or equal to self.annotSpecificityFilter, or if filter is OFF
            if self.annotSpecificityFilter == -1 and self.transcriptSpecificityFilter == -1:
                filteredEnrichmentResults.append( enrich )
            elif self.annotSpecificityFilter == -1 and len( lnc_dict[ transcriptID] ) <= self.transcriptSpecificityFilter:
                filteredEnrichmentResults.append( enrich )
                
            elif len( annot_dict[annotID] ) <= self.annotSpecificityFilter and self.transcriptSpecificityFilter == -1:
                filteredEnrichmentResults.append( enrich )

            elif len( annot_dict[annotID] ) <= self.annotSpecificityFilter and len( lnc_dict[ transcriptID] ) <= self.transcriptSpecificityFilter:
                filteredEnrichmentResults.append( enrich )
            else:
                notPassingFilters += 1
     
        Logger.get_instance().info( "filter_by_specificity : Filtered out %s enrichments" % ( notPassingFilters ) )
         

        return filteredEnrichmentResults


    # #
    # Apply a "best" RNA enrichment filter. For each complex, retain only the X% best enrichments, based on number of observed interactions.
    def filter_by_observed_interactions(self, filtered_enrichment_results):
        
        #===============================================================================       
        # Build dictionary with enrichments per complex
        #===============================================================================
        
        complexEnrichments = {} # key -> complexID, value -> dict; key -> observed interactions, value -> enrichment data
        
        for enrichment in filtered_enrichment_results[1:]:
            spl = enrichment.split("\t")
            annotID = spl[1]
            observedInteractions = int( spl[2])

            if annotID not in complexEnrichments:
                complexEnrichments[ annotID] = {}
            
            if observedInteractions not in complexEnrichments[ annotID]:
                complexEnrichments[ annotID][ observedInteractions] = []
                
            complexEnrichments[ annotID][ observedInteractions].append( enrichment)


        #===============================================================================       
        # Apply filter
        #===============================================================================

        # Filter enrichments based on higher observed interactions 
        observedFilteredResults = [] # store enrichments that pass this filter
        # write header
        observedFilteredResults.append( filtered_enrichment_results[ 0])
        
        for annotID in complexEnrichments:
            
            # calculate number of enrichments amounting to wanted proportion in this complex
            enrichmentsInComplex = sum([len(complexEnrichments[ annotID][observed]) for observed in complexEnrichments[ annotID]])
            nWantedEnrichments = enrichmentsInComplex / float( self.topEnrichmentsPerComplex)
            
            # sort complex enrichments per higher observed interactions
            sortedEnrichments = sorted( complexEnrichments[ annotID], reverse = 1)
            i = 0 # observed interactions looper
            wantedComplexEnrichments = []  # store enrichments         
            
            # Approach: pick enrichments until wanted complexes, in case of draw (enrichments with same observed interactions) pick all of them
            while len( wantedComplexEnrichments) < nWantedEnrichments:
                 
                currentEnrichments = complexEnrichments[ annotID][ sortedEnrichments[ i]]

                for enrich in currentEnrichments:
                    wantedComplexEnrichments.append( enrich)

                i+=1

            observedFilteredResults.extend( wantedComplexEnrichments)

#             print annotID, enrichmentsInComplex, nWantedEnrichments, len( wantedComplexEnrichments), currentEnrichments

            # we should have = or > enrichments than the wanted number/proportion or enrichments            
            assert len( wantedComplexEnrichments) >= nWantedEnrichments

        Logger.get_instance().info( "filter_by_observed_interactions : Filtered out %s enrichments" % ( len( filtered_enrichment_results) - len( observedFilteredResults)) )
            
        return complexEnrichments, observedFilteredResults


    # #
    # Function to write enrichment results and associated files after all the required filters
    def write_enrichment_results(self, filtered_enrichment_results, dict_pairs, count_real_enrichments, count_random_enrichments):

        #===============================================================================
        # Parse and write filtered enrichment results
        #===============================================================================

        # File with same file as enrichment_results file but filtered based on enrichment_per_rna and other filters
        outFile1 = open( self.outputFolder + FilterEnrichmentResults.REPORT_FILTERED_RNA_ANNOT_RESULTS, "w")

        setRNAs = set()
        setAnnots = set()

        # write header
        outFile1.write( filtered_enrichment_results[0])

        for enrichment in filtered_enrichment_results[1:]:

            # First check if all the results with this RNA should be filtered out or not
            spl = enrichment.split("\t")

            # RNA was not filtered out:    
            txID = spl[0]
            annotID = spl[1]

            setRNAs.add( txID)
            setAnnots.add( annotID)
    
            outFile1.write( enrichment + "\n")
            
        outFile1.close()

        #===============================================================================
        # Write list of enriched RNAs (scaffolding candidates)
        #===============================================================================
        with open( self.outputFolder + FilterEnrichmentResults.REPORT_LIST_RNA_SIGN_ENRICH, "w") as outFile2:
            for rna in setRNAs:
                outFile2.write( rna + "\n")       
           
        #===============================================================================    
        # Write summary of enrichment results
        #===============================================================================

        # File with summary of enrichment results
        outFile3 = open( self.outputFolder + FilterEnrichmentResults.REPORT_ENRICHMENT_SUMMARY, "w")
        
        outFile3.write("Total enrichments\tRandom enrichments\tFiltered enrichements\tFiltered enriched complexes\tFiltered enriched lncRNAs\n")
    
        outFile3.write("%s\t%s\t%s\t%s\t%s\n" % ( count_real_enrichments, count_random_enrichments, len( filtered_enrichment_results[1:]), len( setAnnots), len( setRNAs)) )
    
        outFile3.close()
    

    def run(self):

        # Approach: filters are applied through the filteredEnrichmentResults variable, which changes depending on the filtering applied

        # Read row annotation file, if present
        if self.rowAnnotationFile != "":
            Timer.get_instance().step( "Read row annotation file..")    
            self.rowAnnotation = self.read_annotation_file( self.rowAnnotationFile)

        # Read col annotation file, if present
        if self.colAnnotationFile != "":
            Timer.get_instance().step( "Read column annotation file..")    
            self.colAnnotation = self.read_annotation_file( self.colAnnotationFile)
        
        Timer.get_instance().step( "Read Enrichment per RNA file..")    
        listRNASignificantEnrich, countRealEnrichments, countRandomEnrichments = self.read_enrichment_per_rna_file( )

        Timer.get_instance().step( "Read Enrichment results file..")    
        dictPairs, filteredEnrichmentResults = self.read_enrichment_results_file( listRNASignificantEnrich, countRealEnrichments, countRandomEnrichments )

        Timer.get_instance().step( "Write specificity ranking..")    
        annotDict, lncDict = self.rank_by_specificity( dictPairs, self.filterWarningValue)

        # Filter by specificity, if wanted
        if self.annotSpecificityFilter != -1 or self.transcriptSpecificityFilter != -1:
            Timer.get_instance().step( "Write enrichment results file filtered by specificity..")    
            filteredEnrichmentResults = self.filter_by_specificity( filteredEnrichmentResults, annotDict, lncDict)
        
        # Filter by top enrichments (observed interactions), if wanted
        if self.topEnrichmentsPerComplex != -1:
            Timer.get_instance().step( "Write enrichment results file filtered by top enrichments per complex..")    
            _, filteredEnrichmentResults = self.filter_by_observed_interactions ( filteredEnrichmentResults)

        Timer.get_instance().step( "Write filtered enrichment results file..")    
        self.write_enrichment_results( filteredEnrichmentResults, dictPairs, countRealEnrichments, countRandomEnrichments)

            
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
        parser.add_argument('enrichmentPerRNAFile', metavar='enrichmentPerRNAFile', type=str,
                             help='EnrichmentAnalysisStrategy output with one line per RNA, comparing observed number of enrichments versus control.')
        parser.add_argument('enrichmentResultsFile', metavar='enrichmentResultsFile', type=str,
                             help='EnrichmentAnalysisStrategy output with one line per RNA-annotation pair, with the results of enrichment test.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Folder where to write all output files.')
        parser.add_argument('--matrixValueColumn', metavar='matrixValueColumn', type=int, default = 7,
                             help='Column in enrichmentResultsFile to use as value for matrix output file.')
        parser.add_argument('--filterWarningColumn', metavar='filterWarningColumn', type=int, default = 1,
                             help='Whether to exclude enrichments from enrichmentResultsFile based on the warning column ( col 5, 0-based).')
        parser.add_argument('--filterWarningValue', metavar='filterWarningValue', type=str, default = "1.0",
                             help='Value to output in matrix when result has warning flag on. E.g. 1.0 if using p-values, 0 if using significant boolean.')        
        parser.add_argument('--minimumRatio', metavar='minimumRatio', type=str, default = "OFF",
                             help='Float with minimum ratio value between real number of enrichments and random number of enrichments, if ratio is below given value, all enrichments of that RNA will be excluded. Default = "OFF".')
        parser.add_argument('--rowAnnotationFile', metavar='rowAnnotationFile', type=str, default = "",
                             help='TSV file with per row of matrix file (default = RNA). Can have several annotations for same transcript, one per line. E.g. transcriptID\tannotation. Will output file with annotation sorted same way as matrix.')
        parser.add_argument('--colAnnotationFile', metavar='colAnnotationFile', type=str, default = "",
                             help='TSV file with per column of matrix file (default = protein group). Can have several annotations for same protein group, one per line. E.g. annotationID\tannotation. Will output file with annotation sorted same way as matrix.')
        parser.add_argument('--maskMultiple', metavar='maskMultiple', type=int, default = 1,
                             help='Whether to mask annotations when having more than one annotation (val = 1), or display all annotations separated by comma (val = 0). (default = 1).')
        parser.add_argument('--noAnnotationTag', metavar='noAnnotationTag', type=str, default = "Other",
                             help='Text to write for the transcripts that are not in provided annotation files. Default = "Other"')
        parser.add_argument('--noAnnotationFilter', metavar='noAnnotationFilter', type=int, default = 0,
                             help='If on, columns and rows with no annotation will be filtered out from the matrix file (and associated files). Default = 0 (OFF)')
        parser.add_argument('--annotSpecificityFilter', metavar='annotSpecificityFilter', type=int, default = -1,
                             help='Filter out enrichments where the annotation has enrichment to more than X transcripts. Default = -1 (OFF)')
        parser.add_argument('--transcriptSpecificityFilter', metavar='transcriptSpecificityFilter', type=int, default = -1,
                             help='Filter out enrichments where the transcript has enrichment to more than X annotations. Default = -1 (OFF)')
        parser.add_argument('--minimumProteinInteraction', metavar='minimumProteinInteraction', type=int, default = -1,
                             help='Minimum number of proteins in a given annotation with positive interactions for enrichment to be considered. Default = -1 (OFF)')
        parser.add_argument('--topEnrichmentsPerComplex', metavar='topEnrichmentsPerComplex', type=int, default = -1,
                             help='Percentage (0-100) of best enrichments to keep for each complex. Best enrichments defined as more number of observed interactions. In case of draw, keep all drawing enrichments. Default = -1 (OFF)')
           
        # Gets the arguments
        args = parser.parse_args( ) 

        #===============================================================================
        # Argument verification
        #===============================================================================

        if args.minimumRatio != "OFF":
            args.minimumRatio = float( args.minimumRatio)

        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        filterEnrichmentResults = FilterEnrichmentResults(args.enrichmentPerRNAFile, args.enrichmentResultsFile, args.outputFolder, args.matrixValueColumn, \
                                                          args.filterWarningColumn, args.filterWarningValue, args.minimumRatio, args.rowAnnotationFile, \
                                                          args.colAnnotationFile, args.maskMultiple, args.noAnnotationTag, args.noAnnotationFilter, \
                                                          args.annotSpecificityFilter, args.transcriptSpecificityFilter, args.minimumProteinInteraction, \
                                                          args.topEnrichmentsPerComplex)
        
        filterEnrichmentResults.run()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

