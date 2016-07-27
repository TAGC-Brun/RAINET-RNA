import sys
import os
import argparse

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data.RNA import RNA

#===============================================================================
# Started 25-July-2016 
# Diogo Ribeiro
# Based on RBPScore.py
DESC_COMMENT = "Script to attribute annotations to RNAs of a processed catrapid file."
SCRIPT_NAME = "LncRNAScore.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read file with RNA annotation
# 2) Read processed catRAPID file, add annotation
#===============================================================================

class LncRNAScore(object):
    
    # database for which we want cross references
#    XREF_SOURCE_DB = "Ensembl_PRO"
    ID_MAPPING_OUTPUT_FILE = "/transcripts_ensembl_id.tsv"
    ANNOTATION_OUTPUT_FILE = "/annotated_interactions.tsv"
    SEVERAL_ANNOTATION_TAG = "Overlapping_annotations"
       
    def __init__(self, annotation_file, catrapid_file, rainet_db, output_folder, mask_multiple, annotation_column, id_column, no_annotation_tag, cross_reference):

        self.annotationFile = annotation_file
        self.catRAPIDFile = catrapid_file
        self.rainetDB = rainet_db
        self.outputFolder = output_folder
        self.maskMultiple = mask_multiple
        self.annotationColumn = annotation_column
        self.idColumn = id_column
        self.noAnnotationTag = no_annotation_tag
        self.crossReference = cross_reference

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)


    # #
    # Use RAINET DB to retrieve RNA cross references
    def rna_cross_references(self):
        
        # # Get all geneIDs
        # produce dictionary where key is gene ID / gene name and value the transcriptIDs
        query = self.sql_session.query( RNA.transcriptID, RNA.geneID, RNA.externalGeneName ).all()
 
        rnaGeneReference = {} # key -> gene ID, val -> set of transcriptIDs      
        for transcriptID, geneID, externalID in query:     
 
            if geneID not in rnaGeneReference:
                rnaGeneReference[ geneID] = set()
            rnaGeneReference[ geneID].add( str( transcriptID))

            if externalID not in rnaGeneReference:
                rnaGeneReference[ externalID] = set()                 
            rnaGeneReference[ externalID].add( str( transcriptID))
             
        self.rnaGeneReference = rnaGeneReference


    # #
    # Read list of annotation per RNA. Can use cross references to map IDs.
    def read_annotation_file( self):

        #=======================================================================
        # Example file
        #
        # ENST00000545914 cytoplasmic
        # ENST00000304751 nuclear
        # DANCR   Ji2016
        # ENSG00000188825.9       Ji2016
        #=======================================================================

        #=======================================================================
        # initialising
        #=======================================================================

        transcriptAnnotation = {} # Key -> transcript ensemblID, value -> set of annotations

        lineCounter = 0

        #=======================================================================
        # read file
        #=======================================================================

        with open( self.annotationFile, "r") as inFile:
            inFile.readline() # skip header

            notFound = set()
            
            for line in inFile:
                line = line.strip()
                lineCounter+=1
                
                spl = line.split( "\t")
                
                txID = spl[ self.idColumn]
                
                # select column to use as annotation                
                annotationItem = spl[ self.annotationColumn]

                txList = []  
                                              
                if self.crossReference:
                    # if user defines his input dataset as needing cross referencing
                    if txID in self.rnaGeneReference:
                        # e.g. XLOC_000019
                        txList = self.rnaGeneReference[ txID]
                    elif txID.startswith( "ENSG") and "." in txID: 
                        # process gene ID to exclude the .1 etc #e.g. ENSG00000267565.1
                        txID = txID.split( ".")[0]
                        if txID in self.rnaGeneReference:
                            txList = self.rnaGeneReference[ txID]
                    elif txID.startswith( "ENST"):
                        # if actually a normal ENST ID, no need to process
                        txList.append( txID)
                    else:
                        # cannot map ID
                        notFound.add( txID)
                        continue
                else:
                    # a normal ENST ID
                    if not txID.startswith("ENST"):
                        raise RainetException("read_annotation_file: transcript ID is not ENST*:", txID)
                    txList.append( txID)

                for txID in txList:
                    # there can be several txID if using cross references
                    if txID not in transcriptAnnotation:
                        transcriptAnnotation[ txID] = set()
                    transcriptAnnotation[ txID].add( annotationItem)

            print "read_annotation_file: number of entries read:", lineCounter
            print "read_annotation_file: number of transcript IDs not found: ", len( notFound)
            print "read_annotation_file: number of transcripts with annotation:", len( transcriptAnnotation)

        return transcriptAnnotation


    # #
    # Read processed catRAPID file and write output, with added annotation
    def read_catrapid_file( self, transcript_annotation):

        #=======================================================================
        # Example file
        # F8VUJ3  10.7053463052
        # P20933  11.3934851506
        #=======================================================================
  
        #=======================================================================
        # initialising
        #=======================================================================
  
        lineCount = 0
   
        outFile = open( self.outputFolder + LncRNAScore.ANNOTATION_OUTPUT_FILE, "w")
          
        #=======================================================================
        # read file
        #=======================================================================

        with open( self.catRAPIDFile, "r") as inFile:

            # write same header as input file, adding "annotation"
            outFile.write( inFile.readline().strip() + "\tannotation\n")

            for line in inFile:
                line = line.strip()
                lineCount += 1
  
                spl = line.split("\t")
                  
                txID = spl[ 0]

                outputText = ""
                
                if txID in transcript_annotation:
                    # if there is several annotations and we only want information for proteins with a single domain
                    if len( transcript_annotation[ txID]) > 1 and self.maskMultiple == 1:
                        annotation = LncRNAScore.SEVERAL_ANNOTATION_TAG
                    # if wanting information for all annotations
                    elif len( transcript_annotation[ txID]) > 0:
                        # get all annotations
                        annotation = ",".join( list(transcript_annotation[ txID]) )
                    else:
                        raise RainetException("read_catrapid_file: Annotation information is incorrect. ", transcript_annotation[ txID])
                # if there is no annotation for protein
                else:
                    annotation = self.noAnnotationTag

                outputText = "%s\t%s\n" % (line, annotation)

                outFile.write( outputText)
    
        outFile.close()
  
        print "read_catrapid_file: read %s lines.." % lineCount


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
                             help='TSV file with annotation per transcript. No header. Can have several annotations for same transcript, one per line. E.g. transcriptID\tannotation.')
        parser.add_argument('catRAPIDFile', metavar='catRAPIDFile', type=str,
                             help='rnaInteractions.tsv output file from ReadCatrapid.py. E.g. transcriptID\tmean_score')
        parser.add_argument('rainetDB', metavar='rainetDB', type=str, help='Path to RAINET database to be used.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
        parser.add_argument('--maskMultiple', metavar='maskMultiple', type=int, default = 1,
                             help='Whether to mask rna annotations when having more than one annotation (val = 1), or display all annotations separated by comma (val = 0). (default = 1).')
        parser.add_argument('--annotationColumn', metavar='annotationColumn', type=int, default = 1,
                             help='Which column in the input annotation file to process. 0-based.')
        parser.add_argument('--idColumn', metavar='idColumn', type=int, default = 0,
                             help='Which column in the input annotation file to processed as transcript ID. 0-based.')
        parser.add_argument('--noAnnotationTag', metavar='noAnnotationTag', type=str, default = "Other",
                             help='Text to write for the transcripts that are not in provided annotation files. Default = "Other"')
        parser.add_argument('--crossReference', metavar='crossReference', type=str, default = 0,
                             help='Whether to use cross references for transcript ID mapping or not.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        # init
        run = LncRNAScore( args.annotationFile, args.catRAPIDFile, args.rainetDB, args.outputFolder, 
                              args.maskMultiple, args.annotationColumn, args.idColumn, args.noAnnotationTag, args.crossReference)

        # build cross references
        if run.crossReference:
            Timer.get_instance().step( "Reading cross references..")    
            run.rna_cross_references()
        
        # read annotations file
        Timer.get_instance().step( "Reading annotation file..")    
        rnaAnnotation = run.read_annotation_file( )      

        # read catrapid file and write output
        Timer.get_instance().step( "Reading catrapid interaction file..")    
        run.read_catrapid_file( rnaAnnotation)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

