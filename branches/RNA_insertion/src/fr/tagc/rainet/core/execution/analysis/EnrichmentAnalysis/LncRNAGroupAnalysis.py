import sys
import os
import argparse
import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

#===============================================================================
# Started 28-Dec-2016 
# Diogo Ribeiro
# Based on LncRNAScore.py
DESC_COMMENT = "Script to map attributes from Mukherjee2016 to groups of lncRNAs."
SCRIPT_NAME = "LncRNAGroupAnalysis.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read file with RNA annotation, gene IDs
# 2) Read Mukherjee2016 file with data for lncRNAs
# 3) Output into R-readable format
#===============================================================================

#===============================================================================
# Processing notes:
# 1) A category including all RNAs in Mukherjee2016 is created while reading its file
#===============================================================================


class LncRNAGroupAnalysis(object):

    #=======================================================================
    # Constants
    #=======================================================================
    
    ANNOTATION_FILE_ID_COLUMN = 0
    ANNOTATION_FILE_ANNOTATION_COLUMN = 1

    DATA_FILE_ID_COLUMN = 0
    DATA_FILE_ANNOTATION_COLUMN = 9

    OUTPUT_FILE = "lncRNA_group_analysis.tsv"
    
    ALL_MRNA_ANNOTATION = "0-All_mRNAs"
    ALL_LNCRNA_ANNOTATION = "1-All_lncRNAs"


           
    def __init__(self, annotationFile, dataFile, outputFolder, dataColumns, useMRNA):

        self.annotationFile = annotationFile
        self.dataFile = dataFile
        self.outputFolder = outputFolder
        try:
            self.dataColumns = [] 
            sp = dataColumns.split(",")
            for s in sp:
                self.dataColumns.append( int( s))
        except:
            raise RainetException("LncRNAGroupAnalysis.__init__: data column input in wrong format:", dataColumns)
        self.useMRNA = useMRNA

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)



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
        # The same gene can have several annotations

        #=======================================================================
        # initialising
        #=======================================================================

        transcriptAnnotation = {} # Key -> transcript ensemblID, value -> set of annotations
        groupTranscripts = {} # Key -> annotation, value -> set of transcripts

        lineCounter = 0

        #=======================================================================
        # read file
        #=======================================================================

        with open( self.annotationFile, "r") as inFile:
            
            for line in inFile:
                line = line.strip()
                lineCounter+=1
                
                spl = line.split( "\t")
                
                geneID = spl[ LncRNAGroupAnalysis.ANNOTATION_FILE_ID_COLUMN]
                
                # select column to use as annotation                
                annotationItem = spl[ LncRNAGroupAnalysis.ANNOTATION_FILE_ANNOTATION_COLUMN]

                if not geneID.startswith( "ENSG"):
                    raise RainetException("read_annotation_file: entry is not ENSG*:", geneID)

                if "." in geneID:
                    geneID = geneID.split( ".")[0]

                # storing tx as key
                if geneID not in transcriptAnnotation:
                    transcriptAnnotation[ geneID] = set()
                transcriptAnnotation[ geneID].add( annotationItem)
                
                # storing annotation as key
                if annotationItem not in groupTranscripts:
                    groupTranscripts[ annotationItem] = set()
                groupTranscripts[ annotationItem].add( geneID)

            print "read_annotation_file: number of entries read:", lineCounter
            print "read_annotation_file: number of transcripts with annotation:", len( transcriptAnnotation)
            print "read_annotation_file: number of annotations:", len( groupTranscripts)

        self.transcriptAnnotation = transcriptAnnotation
        self.groupTranscripts = groupTranscripts
        
        for group in sorted(groupTranscripts):
            print group, len( groupTranscripts[ group])


    # #
    # Read Mukherjee 2016 file with data
    def read_data_file(self):

        #=======================================================================
        # Example file
        #
        # Gene    Syn     Proc    Deg     CytNuc  PolyCyt TrP     Copies  Exon    Annotation      Cluster Host    Complex
        # ENSG00000005206.12      0.3240500888    -0.0260844809   0.1373502068    -0.5552417614   -0.2815917912   0.6640126412    0.2623901975    MultiExon       lncRNA  c3      None    processed_transcript
        # ENSG00000006062.9       0.1118696177    -0.0129556703   0.3003516672    -0.4050632081   0.0920502949    -0.5617828392   -0.0963797176   MultiExon       lncRNA  c4      None    processed_transcript
        # ENSG00000031544.10      -1.050910308    -0.254916842    0.9567499553    -0.9364242934   -0.1898011997   -2.9665750821   -1.9313555304   MultiExon       lncRNA  c7      None    processed_transcript
        #=======================================================================
        # Note the gene ID has a value after the "."
        # Some classifications are as floats others as strings.
        

        #=======================================================================
        # Output file, a melted file
        #
        # Gene    Group  Metric Value
        # ENSG00000005206 Predicted Syn 0.3240500888
        # ENSG00000005206 Predicted Proc -0.0260844809
        # ENSG00000006062 Interacting Syb 0.1118696177

        outFile = open( self.outputFolder + "/" + LncRNAGroupAnalysis.OUTPUT_FILE, "w")

        # write header
        outFile.write("Gene\tGroup\tMetric\tValue\n")

        numbersPerGroup = {} # key -> group, value -> count of transcripts
        numbersPerGroup[ LncRNAGroupAnalysis.ALL_LNCRNA_ANNOTATION] = 0
        if self.useMRNA:
            numbersPerGroup[ LncRNAGroupAnalysis.ALL_MRNA_ANNOTATION] = 0           

        #=======================================================================
        # read input file and write output
        #=======================================================================

        table = pd.read_table( self.dataFile, header = 0, sep = "\t", skip_blank_lines = True)
        
        columnNames = list(table.columns.values)

        newTable = table[ :]
        
        for index, gene in newTable.iterrows():
            geneID = gene[LncRNAGroupAnalysis.DATA_FILE_ID_COLUMN]

            ## process geneID
            if not geneID.startswith( "ENSG"):
                raise RainetException("read_data_file: entry is not ENSG*:", geneID)

            # Note: some entries contain ENSGR*, this is a small modification due to chromosome Y/X, it can safely be changed to ENSG0*
            if geneID.startswith( "ENSGR"):
                geneID = geneID.replace("ENSGR","ENSG0")

            if "." in geneID:
                geneID = geneID.split( ".")[0]


            # if gene has annotation
            if geneID in self.transcriptAnnotation:
                # for each of its annotations, write a line
                for annotation in self.transcriptAnnotation[ geneID]:
                    
                    for metric in self.dataColumns:
                        outFile.write( "%s\t%s\t%s\t%s\n" % (geneID, annotation, columnNames[ metric], gene[ metric]))

                    if annotation not in numbersPerGroup:
                        numbersPerGroup[ annotation] = 0
                    numbersPerGroup[ annotation]+= 1

            # if mRNA
            if gene[ LncRNAGroupAnalysis.DATA_FILE_ANNOTATION_COLUMN] == "protein_coding":
                if self.useMRNA:
                    # add to mRNA category
                    numbersPerGroup[ LncRNAGroupAnalysis.ALL_MRNA_ANNOTATION]+= 1
                    for metric in self.dataColumns:
                        outFile.write( "%s\t%s\t%s\t%s\n" % (geneID, LncRNAGroupAnalysis.ALL_MRNA_ANNOTATION, columnNames[ metric], gene[ metric]))            
            elif gene[ LncRNAGroupAnalysis.DATA_FILE_ANNOTATION_COLUMN] == "lncRNA":
                # add lncRNA to all lncRNA group regardless of its existence in our annotations
                numbersPerGroup[ LncRNAGroupAnalysis.ALL_LNCRNA_ANNOTATION]+= 1
                for metric in self.dataColumns:
                    outFile.write( "%s\t%s\t%s\t%s\n" % (geneID, LncRNAGroupAnalysis.ALL_LNCRNA_ANNOTATION, columnNames[ metric], gene[ metric]))
            else:
                # neither lncRNA nor mRNA
                continue

        outFile.close()

        print "read_data_file: number of lines in input data:", len(newTable)
        print "read_data_file: number of lncRNAs per group", numbersPerGroup



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
                             help='TSV file with annotation per transcript (gene). No header. Can have several annotations for same transcript, one per line. E.g. transcriptID\tannotation.')
        parser.add_argument('dataFile', metavar='dataFile', type=str,
                             help='File with data per lncRNA from Mukherjee2016. Header is important. Already filtered for lncRNAs.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
        parser.add_argument('--dataColumns', metavar='dataColumns', type=str, default = "1,2,3,4,5,7,10",
                             help='Which 0-based columns in the input data file we want to process. At least the gene ID column needs to be included and as the first in list. Give attribute as comma-separated.')
        parser.add_argument('--useMRNA', metavar='useMRNA', type=int, default = 1,
                             help='Whether to create protein_coding category, if available on file.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        # init
        run = LncRNAGroupAnalysis( args.annotationFile, args.dataFile, args.outputFolder, args.dataColumns, args.useMRNA)
       
        # read annotations file
        Timer.get_instance().step( "Reading annotation file..") 
        run.read_annotation_file( )      

        # read data file and write output
        Timer.get_instance().step( "Reading data file..")    
        run.read_data_file()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

