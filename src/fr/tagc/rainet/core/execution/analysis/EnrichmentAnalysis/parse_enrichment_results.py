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
SCRIPT_NAME = "parse_enrichment_results.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read enrichment_per_rna file, get list of RNAs to filter out based on input criteria.
# 2) Rewrite enrichment_results file based on these filters, output also enrichment matrix file.
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Only writing to file annotations with at least one enrichment
# 2) 
#===============================================================================


#===============================================================================
# File constants
#===============================================================================

REPORT_LIST_RNA_SIGN_ENRICH = "list_RNA_above_random.txt"
REPORT_FILTERED_RNA_ANNOT_RESULTS = "enrichment_results_filtered.tsv"
REPORT_RNA_ANNOT_RESULTS_MATRIX = "enrichment_results_filtered_matrix.tsv"
REPORT_MATRIX_ROW_ANNOTATION = "matrix_row_annotation.tsv"
REPORT_MATRIX_COL_ANNOTATION = "matrix_col_annotation.tsv"

SEVERAL_ANNOTATION_TAG = "Overlapping_annotations"

# #
# Read RNA enrichment file, get list of RNAs with significantly more enrichments compared to control.
def read_enrichment_per_rna_file( enrichment_per_rna_file, minimum_ratio):
       
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

    with open( enrichment_per_rna_file, "r") as inFile:
        
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

            if minimum_ratio != "OFF":
                # filter also by ratio between real and random
                randomSign = float( randomSign)
                
                if randomSign != 0: # to avoid division by Zero
                    ratio = float( observedSign) / randomSign
                else:
                    ratio = 0

                if ratio < minimum_ratio:
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
    Logger.get_instance().info( "read_enrichment_per_rna_file : Total number of RNAs: %s " % len( poolRNA) )
    Logger.get_instance().info( "read_enrichment_per_rna_file : (before filtering) Observed mean: %.1f Random mean: %.1f Mean difference: %.1f" % ( np.mean( observedValues), np.mean( randomValues), np.mean( diffValues)) )

    Logger.get_instance().info( "read_enrichment_per_rna_file : Number of RNAs with enrichment significantly above random: %s " % countAboveRandom )
    if minimum_ratio != "OFF":
        Logger.get_instance().info( "read_enrichment_per_rna_file : Number of RNAs passing real/random ratio: %s " % countRatioPassed )
    Logger.get_instance().info( "read_enrichment_per_rna_file : Number of RNAs passing filters: %s " % len( listRNASignificantEnrich) )
    Logger.get_instance().info( "read_enrichment_per_rna_file : (after filtering) Observed mean: %.1f Random mean: %.1f" % ( np.mean( observedValuesAfter), np.mean( randomValuesAfter)) )

    # average number of tests per RNA, versus control

    # Write list of enriched RNAs to file
    with open( REPORT_LIST_RNA_SIGN_ENRICH, "w") as outFile:
        for rna in listRNASignificantEnrich:
            outFile.write( rna + "\n")       

    return listRNASignificantEnrich

    
# #
# Read RNA-annotation enrichment file (results), filter out results for RNAs that are not significantly enriched over random control, produce output files.
def read_enrichment_results_file(enrichment_results_file, list_rna_significant_enrich, matrix_value_column, filter_warning_column, filter_warning_value,
                                 row_annotation, col_annotation, mask_multiple, no_annotation_tag):
    
    # Example format:
    # transcriptID    annotID number_observed_interactions    number_possible_interactions    total_interacting_proteins      warning pval    corrected_pval  sign_corrected
    # ENST00000230113 344     12      15      7515    0       5.0e-01 7.3e-01 0   

    #===============================================================================
    # Output files
    #===============================================================================
        
    # File with same file as enrichment_results file but filtered based on enrichment_per_rna
    outFile1 = open( REPORT_FILTERED_RNA_ANNOT_RESULTS, "w")

    # File with enrichment_results file data in matrix format.
    # only for RNAs that pass enrichment_per_rna filter
    outFile2 = open( REPORT_RNA_ANNOT_RESULTS_MATRIX, "w")

    #===============================================================================
    # Read file, apply filter, write output
    #===============================================================================

    # initialise some variables
    excludedByRNA = 0
    setRNAs = set()
    setAnnots = set()
    dictPairs = {}
    annotValues = {} # key -> annot ID, val -> list of values

    with open( enrichment_results_file, "r") as inFile:
        
        outFile1.write( inFile.readline()) # transport header

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
            warningFlag = int( spl[5])
            signFlag = int( spl[8])
            
            value = spl[ matrix_value_column] # can be p-value, corrected p-value or significant boolean

            # if filtering is on
            if filter_warning_column:
                # If value is determined not significant OR is tagged with a warning (for several reasons), fill it with constant value. 
                if warningFlag or signFlag == 0:
                    value = filter_warning_value
                else:
                    outFile1.write( line + "\n")

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
        if annotValues[ annotID].count( filter_warning_value) == len( annotValues[ annotID]):
            nonEnrichedAnnotations.add( annotID)

    Logger.get_instance().info( "read_enrichment_results_file : Number of lines filtered out because of enrichment_rna filter: %s" % ( excludedByRNA) )

    assert len( dictPairs) ==  len( setRNAs) * len( setAnnots), "number of pairs should equal number of RNAs times number of annotations"

    outFile1.close()

    #===============================================================================
    # Writing matrix file
    #===============================================================================

    sortedSetRNAs = sorted( setRNAs)   
    sortedSetAnnots = sorted( setAnnots - nonEnrichedAnnotations)

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

    #===============================================================================
    # Writing row annotation file (if applicable)
    #===============================================================================

    if len( row_annotation) > 0:
        # if rowAnnotation option is on
        annotations = [] # contains list of annotations in order
        for rna in sortedSetRNAs:
            if rna in row_annotation:
                # if there is several annotations and we only want information for proteins with a single domain
                if len( row_annotation[ rna]) > 1 and mask_multiple == 1:
                    annotation = SEVERAL_ANNOTATION_TAG
                # if wanting information for all annotations
                elif len(row_annotation[ rna]) > 0:
                    # get all annotations
                    annotation = ",".join( list( row_annotation[ rna]) )
                else:
                    raise RainetException("read_enrichment_results_file: Annotation information is incorrect. ", row_annotation[ rna])
            # if there is no annotation information for current item
            else:
                annotation = no_annotation_tag
        
            annotations.append( annotation)

        # Write file with annotation of rows (default = RNA) based on same sorting as output matrix (for R script)
        outFile3 = open( REPORT_MATRIX_ROW_ANNOTATION, "w")        
        outFile3.write( "\t".join( annotations) + "\n")
        outFile3.close()

    #===============================================================================
    # Writing col annotation file (if applicable)
    #===============================================================================

    if len( col_annotation) > 0:
        # if colAnnotation option is on
        annotations = [] # contains list of annotations in order
        for annotID in sortedSetAnnots:
            if annotID in col_annotation:
                # if there is several annotations and we only want information for proteins with a single domain
                if len( col_annotation[ annotID]) > 1 and mask_multiple == 1:
                    annotation = SEVERAL_ANNOTATION_TAG
                # if wanting information for all annotations
                elif len(col_annotation[ annotID]) > 0:
                    # get all annotations
                    annotation = ",".join( list( col_annotation[ annotID]) )
                else:
                    raise RainetException("read_enrichment_results_file: Annotation information is incorrect. ", col_annotation[ annotID])
            # if there is no annotation information for current item
            else:
                annotation = no_annotation_tag
        
            annotations.append( annotation)

        # Write file with annotation of rows (default = RNA) based on same sorting as output matrix (for R script)
        outFile4 = open( REPORT_MATRIX_COL_ANNOTATION, "w")        
        outFile4.write( "\t".join( annotations) + "\n")
        outFile4.close()
       


# #
# Reads annotation file and returns dictionary of 
def read_annotation_file( annotation_file):
    
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

        # Read row annotation file, if present
        if args.rowAnnotationFile != "":
            Timer.get_instance().step( "Read row annotation file..")    
            rowAnnotation = read_annotation_file( args.rowAnnotationFile)
        else:
            rowAnnotation = {}
    
        # Read col annotation file, if present
        if args.colAnnotationFile != "":
            Timer.get_instance().step( "Read column annotation file..")    
            colAnnotation = read_annotation_file( args.colAnnotationFile)
        else:
            colAnnotation = {}
        
        Timer.get_instance().step( "Read Enrichment per RNA file..")    
        listRNASignificantEnrich = read_enrichment_per_rna_file( args.enrichmentPerRNAFile, args.minimumRatio)

        Timer.get_instance().step( "Read Enrichment results file..")    
        read_enrichment_results_file( args.enrichmentResultsFile, listRNASignificantEnrich, args.matrixValueColumn, args.filterWarningColumn, args.filterWarningValue,
                                      rowAnnotation, colAnnotation, args.maskMultiple, args.noAnnotationTag )

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

