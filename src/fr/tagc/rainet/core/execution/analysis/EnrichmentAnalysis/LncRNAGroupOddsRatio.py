import sys
import os
import argparse
import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.RNA import RNA


#===============================================================================
# Started 16-Jan-2017 
# Diogo Ribeiro
# Based on LncRNAGroupAnalysis.py
DESC_COMMENT = "Script to calculate odds ratio for the overlap between groups of lncRNAs."
SCRIPT_NAME = "LncRNAGroupOddsRatio.py"
#===============================================================================

#===============================================================================
# General plan:
# 0) (option) Read RAINET DB for transcript-gene mapping
# 1) Read lists of lncRNA groups from our in-house analysis
# 2) Read lists of lncRNAs from external sources
# 3) Read list of background lncRNAs (for the statistics)
# 4) Intersect both input lists against background
# 5) Intersect the previous two intersections
# 6) Output odds ratio and pvalues 
#===============================================================================

#===============================================================================
# Processing notes:
# 1) 
#===============================================================================


class LncRNAGroupOddsRatio(object):

    #=======================================================================
    # Constants
    #=======================================================================
    
    ARGUMENT_RAINET_DB_DEFAULT = ""

    ANNOTATION_FILE_ID_COLUMN = 0
    ANNOTATION_FILE_ANNOTATION_COLUMN = 1

           
    def __init__(self, annotationFile, externalFiles, backgroundList, outputFolder, useGenes, rainetDB):

        self.annotationFile = annotationFile
        self.externalFiles = externalFiles
        self.backgroundList = backgroundList
        self.outputFolder = outputFolder
        self.useGenes = useGenes
        self.rainetDB = rainetDB
        
        # check for argument problems, if ID mapping needs to be used, read RAINET DB
        if useGenes:
            if rainetDB == LncRNAGroupOddsRatio.ARGUMENT_RAINET_DB_DEFAULT:
                raise RainetException("LncRNAGroupOddsRatio.__init__: --rainetDB input needs to be provided for converting transcript IDs to gene IDs.")
            else:
                # Build a SQL session to DB
                SQLManager.get_instance().set_DBpath(self.rainetDB)
                self.sql_session = SQLManager.get_instance().get_session()

                self.read_rainet_db()
                
        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)


    # #
    # Get correspondence between transcript IDs and gene IDs using rainetDB
    def read_rainet_db(self):

        rnaCrossReference = {} # key -> transcript ID, val -> ensembl Gene ID
 
        #===============================================================================
        # Query RNA table, which contains gene ID
        #===============================================================================
        query = self.sql_session.query( RNA.transcriptID, RNA.geneID ).all()
        # Note: an gene name points to several ensembl IDs        
         
        # Correspondence should be many-to-1 (a transcript can only have an associated gene, a gene can have many transcripts) 

        geneSet = set()
        
        for transcriptID, geneID in query:
            if transcriptID in rnaCrossReference:
                raise RainetException("LncRNAGroupOddsRatio.read_rainet_db: duplicate transcript ID ")
            
            rnaCrossReference[ transcriptID] = str( geneID)
            
            geneSet.add( geneID)

        Logger.get_instance().info( "LncRNAGroupOddsRatio.read_rainet_db : Number transcripts read %s" % ( len( rnaCrossReference)) )
        Logger.get_instance().info( "LncRNAGroupOddsRatio.read_rainet_db : Number genes read %s" % ( len( geneSet)) )
     
        self.rnaCrossReference = rnaCrossReference
       

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
                
                ID = spl[ LncRNAGroupOddsRatio.ANNOTATION_FILE_ID_COLUMN]
                
                # select column to use as annotation                
                annotationItem = spl[ LncRNAGroupOddsRatio.ANNOTATION_FILE_ANNOTATION_COLUMN]

                # map ID to gene, if needed
                ID = self._transcript_to_gene_mapping( ID)
                if ID == "":
                    continue
                    
                # storing tx as key
                if ID not in transcriptAnnotation:
                    transcriptAnnotation[ ID] = set()
                transcriptAnnotation[ ID].add( annotationItem)
                 
                # storing annotation as key
                if annotationItem not in groupTranscripts:
                    groupTranscripts[ annotationItem] = set()
                groupTranscripts[ annotationItem].add( ID)
 
            Logger.get_instance().info("LncRNAGroupOddsRatio.read_annotation_file: number of entries read: %s" % lineCounter )
            Logger.get_instance().info("LncRNAGroupOddsRatio.read_annotation_file: number of transcripts with annotation: %s" % len( transcriptAnnotation) )
            Logger.get_instance().info("LncRNAGroupOddsRatio.read_annotation_file: number of annotations: %s" % len( groupTranscripts) )
         
        for group in sorted(groupTranscripts):
            Logger.get_instance().info("LncRNAGroupOddsRatio.read_annotation_file %s: %s" % (group, len( groupTranscripts[ group]) ) )
    
 
        self.transcriptAnnotation = transcriptAnnotation
        self.groupTranscripts = groupTranscripts


    # #
    # Read each of the provided external files
    def read_external_files(self):
        
        
        # Dictionary which will contain lists of transcripts
        externalLists = {} # key -> name of list (from file name), value -> list of transcripts
        
        # Read file containing list of files
    
        with open( self.externalFiles, "r") as inFile:
            for line in inFile:
                line = line.strip()
                
                # parse name of list
                fileName = line.split("/")[-1]
                fileName = fileName.split( ".")[0]
                
                if fileName not in externalLists:
                    externalLists[ fileName] = set()               
                
                # read the actual file with list of transcripts
                with open( line, "r") as listFile:
                    for ID in listFile:
                        ID = ID.strip()

                        # map ID to gene, if needed
                        ID = self._transcript_to_gene_mapping( ID)
                        if ID == "":
                            continue

                        externalLists[ fileName].add( ID)



        Logger.get_instance().info("LncRNAGroupOddsRatio.read_external_files: number of lists read: %s" % len( externalLists) )
        
        for lst in sorted( externalLists):
            Logger.get_instance().info("LncRNAGroupOddsRatio.read_external_files %s: %s" % (lst, len( externalLists[ lst]) ) )


        self.externalLists = externalLists


    # #
    # Read background list of transcripts, filters out other lists.
    def read_background_list(self):
        
        #=======================================================================
        # Read file with list of transcripts 
        #=======================================================================

        backgroundTranscripts = set() 

        with open( self.backgroundList, "r") as inFile:
            for line in inFile:
                ID = line.strip()

                # map ID to gene, if needed
                ID = self._transcript_to_gene_mapping( ID)
                if ID == "":
                    continue
                
                backgroundTranscripts.add( ID)

        Logger.get_instance().info("LncRNAGroupOddsRatio.read_background_list: number of transcripts in background: %s" % len( backgroundTranscripts) )

        #=======================================================================
        # Filter for overlap with background transcripts
        #=======================================================================
        
        ### annotation transcripts
                
        annotationTranscripts = self.transcriptAnnotation.keys()        
        filteredAnnotationTranscripts = [ tx for tx in annotationTranscripts if tx in backgroundTranscripts]
        
        #Note: the overlap should be 100%, but to keep the code stable and flexible it is good to confirm it
        Logger.get_instance().info("LncRNAGroupOddsRatio.read_background_list: annotation transcripts overlap with background: %s out of %s" % ( len( filteredAnnotationTranscripts), len( annotationTranscripts) ) )
        
        ### External transcript lists

        filteredExternalLists = {} # key -> name of list (from file name), value -> list of transcripts overlapping with background 

        for lst in self.externalLists:
            
            filteredTxs = [ tx for tx in self.externalLists[ lst] if tx in backgroundTranscripts]
            
            filteredExternalLists[ lst] = filteredTxs
            
            Logger.get_instance().info("LncRNAGroupOddsRatio.read_background_list: %s transcripts overlap with background: %s out of %s" % ( lst, len( filteredExternalLists[ lst]), len( self.externalLists[ lst]) ) )
            

        self.backgroundTranscripts = backgroundTranscripts
        self.filteredExternalLists = filteredExternalLists
        self.filteredAnnotationTranscripts = filteredAnnotationTranscripts


    # #
    # Calculate odds ratio between groups of lncRNAs
    def calculate_odds_ratio(self):

        self.filteredExternalLists
        self.filteredAnnotationTranscripts
        self.backgroundTranscripts
        #TODO: continue here. Do fisher exact test
        pass



    # #
    # Map transcript ID to gene ID based on RAINET DB correspondence, if needed.
    def _transcript_to_gene_mapping(self, ID):
        
        if ID.startswith("ENSG") and self.useGenes == 0:
            raise RainetException("LncRNAGroupOddsRatio._transcript_to_gene_mapping: using Gene IDs in input file but --useGenes option not chosen.")
        elif ID.startswith( "ENST") and self.useGenes:
            # convert transcriptID to geneID                            
            if ID in self.rnaCrossReference:
                ID = self.rnaCrossReference[ ID]
            else:
                Logger.get_instance().warning( "LncRNAGroupOddsRatio._transcript_to_gene_mapping : Transcript ID not present for ID mapping: %s" % ( ID) )
                ID = ""
        else:
            if not ID.startswith("ENS"):
                Logger.get_instance().warning( "LncRNAGroupOddsRatio._transcript_to_gene_mapping : ID not identified: %s" % ( ID) )
                ID = ""                             

        return ID

    # #
    # Run functions in correct order
    def run(self):
        
        Timer.get_instance().step( "Reading annotation file..") 
        run.read_annotation_file( )      

        Timer.get_instance().step( "Reading external files..") 
        run.read_external_files()
        
        Timer.get_instance().step( "Reading background file..") 
        run.read_background_list()

#         # read data file and write output
#         Timer.get_instance().step( "Reading data file..")    
#         run.read_data_file()



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
                             help='TSV file with annotation per transcript or gene. No header. Can have several annotations for same transcript, one per line. E.g. transcriptID\tannotation.')
        parser.add_argument('externalFiles', metavar='externalFiles', type=str,
                             help='File with list of input files, each with list of transcripts/genes to be used to cross our data. Full paths, no header.')
        parser.add_argument('backgroundList', metavar='backgroundList', type=str,
                             help='File with list of transcripts to use has background (i.e. tested transcripts). Both annotationFile and externalFiles transcripts will be filtered by this list. No header.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
        parser.add_argument('--useGenes', metavar='useGenes', type=str, default = 0,
                             help='If on, files containing ENST* IDs will be mapped to ENSG* IDs. Note: --rainetDB parameter is mandatory.')
        parser.add_argument('--rainetDB', metavar='rainetDB', type=str, default = LncRNAGroupOddsRatio.ARGUMENT_RAINET_DB_DEFAULT,
                             help='Rainet database to be used for ID mapping.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        # init
        run = LncRNAGroupOddsRatio( args.annotationFile, args.externalFiles, args.backgroundList, args.outputFolder, args.useGenes, args.rainetDB)

        # run functions in order
        run.run()
       
        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

