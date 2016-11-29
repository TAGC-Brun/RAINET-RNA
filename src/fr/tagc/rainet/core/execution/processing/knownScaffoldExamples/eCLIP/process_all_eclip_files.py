import sys
import os
import argparse
import glob

# import numpy as np
import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil

from constants import *

# Needs bedtools v2.17.0

#===============================================================================
# Started 24-Nov-2016 
# Diogo Ribeiro
DESC_COMMENT = "Wrapper script for ENCODE eCLIP dataset. To filter, merge replicates, map against gene models, produce interactions file."
SCRIPT_NAME = "process_all_eclip_files.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read metadata file, have group replicates, match identifiers
# 2) Read bed files, apply filtering
# 3) Merge replicate bed files
# 4) Map peaks to transcripts
# 5) Produce interaction file
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Filtering of peaks based on values of columns is done using another script
# 2) Merge of replicates is done by bedtools intersect default behaviour. i.e. keeping only the coordinates that overlap between the two files
# 3) Mapping peaks to transcript is done by bedtools intersect with -wa flag, i.e. keeping the entry corresponding to the transcript models that overlap (even if by only one base) with a peak
#===============================================================================


FILTER_ECLIP_FILES_SCRIPT_PATH = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/knownScaffoldExamples/eCLIP/filter_eclip_files.py"

TEMP_SUFFIX = "_temp"
COMBINE_SUFFIX = "_intersect.bed"
MAP_SUFFIX = "_transcripts.bed"
INTERACTION_FILE_NAME = "eclip_interactions.out"

# On metadata file, which assembly to use (others will be filtered out)
WANTED_ASSEMBLY = "GRCh38"

# Assert that all samples have exactly two replicates
NUMBER_REPLICATES = 2

# List of cell lines names
CELL_LINES = ["K562", "HepG2"]

# ID mapping columns
ORIGINAL_ID = 0
NEW_ID = 2

# #
# Read, filter and write filtered bed file
def read_bed_folder( bedFolder, outputFolder, minPval, minFC):

    filesToProcess = glob.glob( bedFolder + "/*.bed" )

    print "read_bed_folder: Processing %s files.." % len( filesToProcess)
    
    newFilesToProcess = [] # new list of files to process after this step

    count = 0
    for fi in filesToProcess:
        count += 1
        if count % 50 == 0:
            print "read_bed_folder: Processed %s files.." % count
        
        fileName = fi.split("/")[-1]

        # create backup of original file into output folder
        cmd = "cp %s %s/%s%s" % ( fi, outputFolder, fileName, TEMP_SUFFIX)
        SubprocessUtil.run_command( cmd, verbose = 0)

        # new filename
        newFi = "%s/%s%s" % ( outputFolder, fileName, TEMP_SUFFIX)

        # run command to filter file. Output will be created in the same folder 
        #workon rainet;     
        cmd = "python %s %s --minPval %s --minFC %s" % ( FILTER_ECLIP_FILES_SCRIPT_PATH, newFi, minPval, minFC)
        SubprocessUtil.run_command( cmd, verbose = 0)

        # keep name of filtered files, we want to process further
        newFilesToProcess.append( newFi + FITLERED_OUTPUT_SUFFIX)

    # remove temp files
    cmd = "rm %s/*_temp" % ( outputFolder)
    SubprocessUtil.run_command( cmd, verbose = 0)

    return newFilesToProcess


# #
# Read metadata file 
def read_metadata_file( metadataFile):
  
    # Example format
    # File accession  File format     Output type     Experiment accession    Assay   Biosample term id       Biosample term name     Biosample type  Biosample life stage    Biosample sex   Biosample organis  Biosample treatments    Biosample subcellular fraction term name        Biosample phase Biosample synchronization stage Experiment target       Antibody accession      Library made from       Library depleted in     Library extraction method       Library lysis method    Library crosslinking method     Experiment date released        Project RBNS protein concentration      Library fragmentation method    Library size range      Biosample Age   Biological replicate(s) Technical replicate     Read length     Mapped read length      Run type        Paired end      Paired with     Derived from    Size    Lab     md5sum  File download URL       Assembly        Platform        Controlled by   File Status     Audit WARNING   Audit INTERNAL_ACTION   Audit NOT_COMPLIANT     Audit ERROR
    # ENCFF066PCT    bed narrowPeak    peaks    ENCSR438GZQ    eCLIP    EFO:0002067    K562    immortalized cell line    adult    female    Homo sapiens                    KHSRP-human    ENCAB531WWS    RNA        see document    see document    ultraviolet irradiation    2015-12-16    ENCODE        see document    175-300    53 year    1    1                        ENCFF405WBU, ENCFF698TGL    360587    Gene Yeo, UCSD    634b9db066421eb062fcf3e7972d8c4b    https://www.encodeproject.org/files/ENCFF066PCT/@@download/ENCFF066PCT.bed.gz    GRCh38            released        experiment not submitted to GEO, NTR assay, mismatched file status        

    # Columns of interest:
    # File accession - 0, file format - 1, experiment accession - 3, Biosample term name - 6, Experiment target - 15, Biological replicate(s) - 28, Assembly - 40
    wantedIndexes = [0,1,3,6,15,28,40]

    # Bed file names are represented as File accession 
    # Assembly, used for filtering wanted assembly mapped against

    proteinDict = {} # key -> protein name (Experiment target), dict. Key -> cell line, value -> all other information (one item for each biological replicate)

    count = 0    
    with open( metadataFile, "r") as inFile:
        header = inFile.readline()
        for line in inFile:
            spl = line.strip().split('\t')
            
            # include only entries using wanted assembly
            
            assembly = spl[ 40]
            if assembly != WANTED_ASSEMBLY:
                continue
            
            proteinName = spl[ 15]
            cellLine = spl[6]

            # confirm format is bed narrowPeak
            assert spl[1] == "bed narrowPeak"
            
            text = []
            for idx in wantedIndexes:
                text.append( spl[ idx] )

            if proteinName not in proteinDict:
                proteinDict[ proteinName] = {}

            if cellLine not in proteinDict[ proteinName]:
                proteinDict[ proteinName][ cellLine] = []
            proteinDict[ proteinName][ cellLine].append( text)

            count += 1

    print "read_metadata_file: Total %s proteins in metadata file" % ( len( proteinDict))
    print "read_metadata_file: Total %s bed files in metadata file (after filter)" % ( count)

    # make data structure with pairs of experiments consisting of the two biological replicates
    pairsOfReplicates = {} # key -> protein+celline tag, val -> pair of biological replicates
    for proteinName in proteinDict:
        for cellLine in proteinDict[ proteinName]:
            assert len( proteinDict[ proteinName][ cellLine]) == NUMBER_REPLICATES 
            tag = proteinName + "_" + cellLine #e.g. CSTF2T-human_K562
            pair = proteinDict[ proteinName][ cellLine][0][0] + "|" + proteinDict[ proteinName][ cellLine][1][0] #e.g. ENCFF457MNN|ENCFF288SMV
            
            pairsOfReplicates[ tag] = pair
    
    assert len( pairsOfReplicates) == count / NUMBER_REPLICATES

    print "read_metadata_file: Total %s pairs of replicates" % len( pairsOfReplicates)

    return proteinDict, pairsOfReplicates


# #
# After bed filtering, use bedtools intersect to combine data from two replicates
def combine_replicates( filesToProcess, pairsOfReplicates, outputFolder):
    
    # create correspondence between File Accession and file path
    fileAccession = {} # key -> File Accession number, val -> file path
    for fi in filesToProcess:
        fileName = fi.split("/")[-1] #e.g. ENCFF138MEW.bed_temp_filtered
        accession = fileName.split(".")[0]
        
        fileAccession[ accession] = fi

    newFilesToProcess = []

    for prot in pairsOfReplicates:
        pair = pairsOfReplicates[ prot]
        #.e.g ENCFF951RJR|ENCFF655RBH
        fileA, fileB = pair.split("|") 
        
        # get file paths
        fileAPath = fileAccession[ fileA]
        fileBPath = fileAccession[ fileB]

        outFile = "%s/%s%s" % (outputFolder, prot, COMBINE_SUFFIX)

        # run bedtools interserct with default parameters
        cmd = "bedtools intersect -a %s -b %s > %s" % ( fileAPath, fileBPath, outFile)
        SubprocessUtil.run_command( cmd, verbose = 0)

        # append processed file name to list
        newFilesToProcess.append( outFile)

    # remove intermediary files
    cmd = "rm %s/*temp%s" % ( outputFolder, FITLERED_OUTPUT_SUFFIX)
    SubprocessUtil.run_command( cmd, verbose = 0)

    print "combine_replicates: Created %s bed intersect files.." % len( newFilesToProcess)
    
    return newFilesToProcess


# #
# After combining data from replicates, use bed intersect against transcript models to map those peaks to human transcripts 
def map_to_transcript( filesToProcess, bedModels, outputFolder):

    newFilesToProcess = []

    for fi in filesToProcess:

        # set output file name
        fileName = fi.split("/")[-1]
        fileName = fileName.replace( COMBINE_SUFFIX, "")
        outFile = "%s/%s%s" % (outputFolder, fileName, MAP_SUFFIX)
        newFilesToProcess.append( outFile)

        # run bedtools intersect        
        # here the -a and -b matter, as we want to keep transcript information in the output file.
        cmd = "bedtools intersect -wa -a %s -b %s > %s" % ( bedModels, fi, outFile)
        SubprocessUtil.run_command( cmd, verbose = 0)

    print "map_to_transcript: Created %s bed mapping files.." % len( newFilesToProcess)

    return newFilesToProcess


# # 
# Read uniprot ID mapping file 
def read_id_mapping( idMappingFile):

    idMapping = {} # key -> proteinName, value -> uniprotAC

    with open( idMappingFile, "r") as inFile:
        for line in inFile:
            spl = line.split( "\t")
            proteinName = spl[ ORIGINAL_ID]
            proteinAC = spl[ NEW_ID]

            if proteinName not in idMapping:
                idMapping[ proteinName] = proteinAC
            else:
                print "read_id_mapping: problem, using two ACs for same protein name %s" % proteinName

    print "read_id_mapping: using ID mapping of %s proteins" % len( idMapping)

    return idMapping

# #
# From peak-transcript mapping, create an interactions file of all the data
def produce_interactions_file( filesToProcess, outputFolder, idMapping):

    # example mapped bed    
    #chr1    149842203       149886641       ENST00000369160.3       .       -       ENSEMBL transcript      .       ID=ENST00000369160.3;Parent=ENSG00000184678.9;gene_id=ENSG00000184678.9;transcript_id=ENST00000369160.3;gene_type=protein_coding;gene_status=KNOWN;gene_name=HIST2H2BE;transcript_type=protein_coding;transcript_status=KNOWN;transcript_name=HIST2H2BE-201;level=3;protein_id=ENSP00000375736.2;transcript_support_level=5;tag=basic,appris_principal_1;havana_gene=OTTHUMG00000012095.1

    # output file:
    # transcriptID|proteinName\tcount_peaks
    
    outFile = open( outputFolder + "/" + INTERACTION_FILE_NAME, "w")
    
    interactionsPerProtein = {} # key -> proteinName, value -> dict. key -> transcriptID, value -> count number of peaks in transcript
    
    # loop files
    for fi in sorted( filesToProcess):
        
        # get protein name from the file name
        proteinName = fi.split("/")[-1].replace( MAP_SUFFIX, "")
        # remove the cell line tag
        for cellLine in CELL_LINES:
            if cellLine in proteinName:
                proteinName = proteinName.replace( cellLine, "").replace("_","")
        # remove the 'human' tag
        if "-human" in proteinName:
            proteinName = proteinName.replace("-human","")

        # if ID mapping was used as option, replace name of protein with its uniprotAC
        if proteinName in idMapping:
            proteinName = idMapping[ proteinName]

        if proteinName not in interactionsPerProtein:
            interactionsPerProtein[ proteinName] = {}
        
        # open and read each bed file to get list of transcripts
        with open( fi,"r") as inFile:
            for line in inFile:
                spl = line.split("\t")
                transcriptID = spl[3]
                
                # remove number after the dot
                if "." in transcriptID:
                    transcriptID = transcriptID.split(".")[0]

                if transcriptID not in interactionsPerProtein[ proteinName]:
                    interactionsPerProtein[ proteinName][ transcriptID] = 0
                
                interactionsPerProtein[ proteinName][ transcriptID] += 1

    # write all interactions
    for prot in sorted( interactionsPerProtein):
        for txID in sorted( interactionsPerProtein[ prot]):
            outFile.write( "%s|%s\t%s\n" % ( txID, prot, interactionsPerProtein[ prot][ txID] ) )

    outFile.close()
    
    # clean output
    os.system( "mkdir %s/processed" % ( outputFolder))
    os.system( "mv %s/*.bed %s/processed" % ( outputFolder, outputFolder))

    print "produce_interactions_file: produced interactions for %s proteins.." % len( interactionsPerProtein)



if __name__ == "__main__":

    try:
    
        # Start chrono
        Timer.get_instance().start_chrono()
        print "STARTING " + SCRIPT_NAME
        
        #===============================================================================
        # Get input arguments
        #===============================================================================
        parser = argparse.ArgumentParser(description= DESC_COMMENT) 
    
        # positional args
        parser.add_argument('bedFolder', metavar='bedFolder', type=str,
                             help='Input folder with all bed files to be processed.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
        parser.add_argument('metadataFile', metavar='metadataFile', type=str,
                             help='Metadata file provided when downloding the bed files.')
        parser.add_argument('bedModels', metavar='bedModels', type=str,
                             help='Bed file with gene models for identifying overlapping transcripts.')
        parser.add_argument('--minPval', metavar='minPval', type=float, default = MIN_PVAL_DEFAULT,
                             help='Peaks below given value will be excluded. Note that provided pvalue is positive log10, the higher the value, the more significant it is. (Default = -1, i.e. "OFF").')
        parser.add_argument('--minFC', metavar='minPval', type=float, default = MIN_FC_DEFAULT,
                             help='Peaks below given value will be excluded. Note that provided fold-change enrichment is positive log2, the higher the value, the higher is the fold change. (Default = -1, i.e. "OFF").')
        parser.add_argument('--idMapping', metavar='idMapping', type=str, default = "",
                             help='If Uniprot ID mapping file provided, final output file with display interactions of protein name in column 0 using the ID in column 2 (0-based).')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        if not os.path.exists( args.outputFolder):
            os.mkdir( args.outputFolder)

        Timer.get_instance().step( "Read and filter bed folder..")            
        filteredFilesToProcess = read_bed_folder( args.bedFolder, args.outputFolder, args.minPval, args.minFC)
        
        Timer.get_instance().step( "Read metadata file..") 
        proteinDict, pairsOfReplicates = read_metadata_file( args.metadataFile)

        Timer.get_instance().step( "Combine biological replicates..") 
        combinedFilesToProcess = combine_replicates( filteredFilesToProcess, pairsOfReplicates, args.outputFolder)

        Timer.get_instance().step( "Map peak to transcript..") 
        mappedFilesToProcess = map_to_transcript( combinedFilesToProcess, args.bedModels, args.outputFolder)

        if args.idMapping != "":
            Timer.get_instance().step( "Read protein ID mapping file..") 
            idMapping = read_id_mapping( args.idMapping)
        else:
            idMapping = {}

        Timer.get_instance().step( "Produce interaction file..") 
        produce_interactions_file( mappedFilesToProcess, args.outputFolder, idMapping)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

