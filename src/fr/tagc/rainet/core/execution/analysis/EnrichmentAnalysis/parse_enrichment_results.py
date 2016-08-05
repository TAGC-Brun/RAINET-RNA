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
DESC_COMMENT = "Script to parse output from EnrichmentAnalysisStrategy"
SCRIPT_NAME = "parse_enrichment_results.py"
#===============================================================================

#===============================================================================
# General plan:
# 1)
# 2)
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Only writing to file annotations with at least one enrichment
# 2) CURRENTLY NOT WORKING IF USING SIGN COLUMN AS INPUT
#===============================================================================


#===============================================================================
# File constants
#===============================================================================

REPORT_LIST_RNA_SIGN_ENRICH = "list_RNA_above_random.txt"
REPORT_FILTERED_RNA_ANNOT_RESULTS = "enrichment_results_significant.tsv"
REPORT_RNA_ANNOT_RESULTS_MATRIX = "enrichment_results_filtered_matrix.tsv"

WARNING_FILTER_VALUE =  1.0 #NA

# #
# Read RNA enrichment file, get list of RNAs with significantly more enrichments compared to control.
def read_enrichment_per_rna_file( enrichment_per_rna_file):
    
    
    # Example format:
    # transcriptID    n_sign_tests_no_warning avg_n_sign_random_no_warning    n_times_above_random    empiricalPvalue significant
    # ENST00000230113 57      50.78   815     1.9e-01 0

    poolRNA = set()

    listRNASignificantEnrich = set()

    observedValues = []
    randomValues = []
    diffValues = []

    with open( enrichment_per_rna_file, "r") as inFile:
        
        inFile.readline() # skip header

        for line in inFile:
            line = line.strip()

            # unpack line            
            txID, observedSign, randomSign, _, _, signFlag = line.split( "\t")
            
            poolRNA.add( txID)

            if signFlag == "1":
                listRNASignificantEnrich.add( txID)
                        
            observedValues.append( float(observedSign))
            randomValues.append( float(randomSign))

            # difference between observed and random number of tests with enrichment
            diffValue = float( observedSign) - float( randomSign)
            diffValues.append( diffValue)

    # Print basic stats
    Logger.get_instance().info( "read_enrichment_per_rna_file : Total number of RNAs: %s " % len( poolRNA) )
    Logger.get_instance().info( "read_enrichment_per_rna_file : Number of RNAs with enrichment significantly above random: %s " % len( listRNASignificantEnrich) )

    Logger.get_instance().info( "read_enrichment_per_rna_file : Observed mean: %.1f Random mean: %.1f Mean difference: %.1f" % ( np.mean( observedValues), np.mean( randomValues), np.mean( diffValues)) )
    # average number of tests per RNA, versus control

    # Write list of enriched RNAs to file
    with open( REPORT_LIST_RNA_SIGN_ENRICH, "w") as outFile:
        for rna in listRNASignificantEnrich:
            outFile.write( rna + "\n")       

    return listRNASignificantEnrich

    
# #
# Read RNA-annotation enrichment file (results), filter out results for RNAs that are not significantly enriched over random control, produce output files.
def read_enrichment_results_file(enrichment_results_file, list_rna_significant_enrich, matrix_value_column, filter_warning_column):
    
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
                    value = WARNING_FILTER_VALUE
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
        if annotValues[ annotID].count( WARNING_FILTER_VALUE) == len( annotValues[ annotID]):
            nonEnrichedAnnotations.add( annotID)

    Logger.get_instance().info( "read_enrichment_results_file : Number of lines filtered out because of enrichment_rna filter: %s" % ( excludedByRNA) )

    assert len( dictPairs) ==  len( setRNAs) * len( setAnnots), "number of pairs should equal number of RNAs times number of annotations"


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
    
    outFile1.close()
    outFile2.close()
    

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
           
        # Gets the arguments
        args = parser.parse_args( ) 

        #===============================================================================
        # Run analysis / processing
        #===============================================================================
    
        Timer.get_instance().step( "Read Enrichment per RNA file..")    
        listRNASignificantEnrich = read_enrichment_per_rna_file( args.enrichmentPerRNAFile)

        Timer.get_instance().step( "Read Enrichment results file..")    
        read_enrichment_results_file( args.enrichmentResultsFile, listRNASignificantEnrich, args.matrixValueColumn, args.filterWarningColumn)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

