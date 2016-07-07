import sys
import os
import argparse

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference

#===============================================================================
# Started 22-May-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to attribute annotations to proteins of a processed catrapid file."
SCRIPT_NAME = "RBPDomainScore.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read file with protein domain or other annotation
# 2) Read processed catRAPID file, add annotation
#===============================================================================

class RBPDomainScore(object):
    
    # database for which we want cross references
    XREF_SOURCE_DB = "Ensembl_PRO"
    UNIPROT_OUTPUT_FILE = "/proteins_uniprotac.tsv"
    ANNOTATION_OUTPUT_FILE = "/annotated_interactions.tsv"
    #NO_ANNOTATION_TAG = "Non-binding"
    SEVERAL_ANNOTATION_TAG = "Overlapping_annotations"
       
    def __init__(self, domain_annotation_file, catrapid_file, rainet_db, output_folder, mask_multiple, annotation_column, id_column, extra_annotation, no_annotation_tag):

        self.annotationFile = domain_annotation_file
        self.catRAPIDFile = catrapid_file
        self.rainetDB = rainet_db
        self.outputFolder = output_folder
        self.maskMultiple = mask_multiple
        self.annotationColumn = annotation_column
        self.idColumn = id_column
        self.extraAnnotation = extra_annotation
        self.noAnnotationTag = no_annotation_tag

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)


    # #
    # Use RAINET DB to retrieve Protein cross references
    def protein_cross_references(self):
        
        # Query all UniProtAC in database
        query = self.sql_session.query( Protein.uniprotAC ).all()        
        uniprotACs = { str(prot[0]) for prot in query}

        # # Get external references
        # produce dictionary where key is xref ID and value the uniprotAC
        # Note: an external ID can point to several uniprot IDs        
        query = self.sql_session.query( ProteinCrossReference.protein_id, ProteinCrossReference.crossReferenceID ).filter( ProteinCrossReference.sourceDB == RBPDomainScore.XREF_SOURCE_DB).all()

        protCrossReference = {} # key -> external ID, val -> set of uniprotACs        
        for uniprotID, externalID in query:     

            if externalID not in protCrossReference:
                protCrossReference[ externalID] = set()
                
            protCrossReference[ externalID].add( str( uniprotID))
            
        return uniprotACs, protCrossReference


    # #
    # Read list of domains per protein. Use cross references
    def read_domain_annotation_file(self, protCrossReference):

        #=======================================================================
        # Example file
        #
        # gene id prot id pfam id domain name     domain class
        # ENSG00000269612 ENSP00000469200 PF00270.26,PF00271.28   DEAD,Helicase_C Classical,Non-classical
        #=======================================================================

        #=======================================================================
        # initialising
        #=======================================================================

        proteinAnnotation = {} # Key -> uniprotAC, value -> set of annotations

        lineCounter = 0

        outFile = open( self.outputFolder + RBPDomainScore.UNIPROT_OUTPUT_FILE, "w")
        # output file header
        outFile.write( "ensemblID\tuniprotacs\tannotations\n")

        #=======================================================================
        # read file
        #=======================================================================

        with open( self.annotationFile, "r") as inFile:
            # skip header
            inFile.readline()

            notFound = set()
            
            for line in inFile:
                line = line.strip()
                lineCounter+=1
                
                spl = line.split( "\t")
                
                pepID = spl[ self.idColumn]
                
                # select column to use as annotation                
                annotationItem = spl[ self.annotationColumn]
                                
                if not pepID.startswith("ENSP"):
                    raise RainetException("read_domain_annotation_file: Peptide ID is incorrect:", pepID)

                # retrieve corresponding uniprotACs
                if pepID in protCrossReference:
                    protIDs = protCrossReference[ pepID]
                else:
                    notFound.add( pepID)
                    continue

                # there can be several domain classes for the same protein
                annotationItems = annotationItem.split( ",")

                outFile.write( "%s\t%s\t%s\n" % ( pepID, ",".join( protIDs), annotationItem ))

                for cla in annotationItems:
                    # there can be several protein uniprotAC for a single Ensembl peptideID
                    for protID in protIDs:
                        if protID not in proteinAnnotation:
                            proteinAnnotation[ protID] = set()
#                         else:
#                             print "read_domain_annotation_file: duplicate annotation line for same protein.", protID, line
                        proteinAnnotation[ protID].add( cla)

            print "read_domain_annotation_file: number of entries read:", lineCounter
            print "read_domain_annotation_file: number of peptide IDs not found: ", len(notFound)
            print "read_domain_annotation_file: number of proteins with annotation:", len(proteinAnnotation)

        outFile.close()

        return proteinAnnotation


    # #
    # Read extra list of annotations per protein uniprotac. Without using cross references.
    # @return proteinAnnotation, with new annotations if case that input file was provided
    # Note: currently not writing into UNIPROT_OUTPUT_FILE.
    def read_extra_annotation_file(self, proteinAnnotation):

        # check if extra annotation file provided
        if self.extraAnnotation != "":
            counter = 0
            countNewProt = 0
            with open( self.extraAnnotation, "r") as inFile:
                for line in inFile:
                    spl = line.strip().split("\t")
                    uniprotac, annot = spl[0], spl[1]
                    
                    # add information to already existing proteins, or create new one                    
                    if uniprotac not in proteinAnnotation:
                        proteinAnnotation[ uniprotac] = set()
                        countNewProt += 1

                    proteinAnnotation[ uniprotac].add( annot)
                    counter += 1
                    
            print "read_extra_annotation_file: number of annotations added with extra file:", counter
            print "read_extra_annotation_file: number of new proteins added with extra file:", countNewProt
                        
            return proteinAnnotation
        
        # if no file provided return object without changes
        else:
            return proteinAnnotation



    # #
    # Read processed catRAPID file and write output, with added annotation
    def read_catrapid_file( self, proteinAnnotation):

        #=======================================================================
        # Example file
        # F8VUJ3  10.7053463052
        # P20933  11.3934851506
        #=======================================================================
  
        #=======================================================================
        # initialising
        #=======================================================================
  
        lineCount = 0
   
        outFile = open( self.outputFolder + RBPDomainScore.ANNOTATION_OUTPUT_FILE, "w")
          
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
                  
                protID = spl[0]
                score = spl[1]

                outputText = ""
                
                if protID in proteinAnnotation:
                    # if there is several annotations and we only want information for proteins with a single domain
                    if len( proteinAnnotation[ protID]) > 1 and self.maskMultiple == 1:
                        annotation = RBPDomainScore.SEVERAL_ANNOTATION_TAG
                    # if wanting information for all annotations
                    elif len( proteinAnnotation[ protID]) > 0:
                        # get all annotations
                        annotation = ",".join( list(proteinAnnotation[ protID]) )
                    else:
                        raise RainetException("read_catrapid_file: Annotation information is incorrect. ", proteinAnnotation[ protID])
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
                             help='TSV file with domain annotation per protein. E.g. gene id prot id pfam id domain name     domain class.')
        parser.add_argument('catRAPIDFile', metavar='catRAPIDFile', type=str,
                             help='proteinInteractions.tsv output file from ReadCatrapid.py. E.g. uniprotac\tmean_score')
        parser.add_argument('rainetDB', metavar='rainetDB', type=str, help='Path to RAINET database to be used.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
        parser.add_argument('--maskMultiple', metavar='maskMultiple', type=int, default = 1,
                             help='Whether to mask protein annotations when having more than one annotation (val = 1), or display all annotations separated by comma (val = 0). (default = 1).')
        parser.add_argument('--annotationColumn', metavar='annotationColumn', type=int, default = 4,
                             help='Which column in the input annotation file to process. 0-based.')
        parser.add_argument('--idColumn', metavar='idColumn', type=int, default = 1,
                             help='Which column in the input annotation file to processed as protein ID. 0-based.')
        parser.add_argument('--extraAnnotation', metavar='extraAnnotation', type=str, default = "",
                             help='File with extended annotation. Column1: uniprotAc, Column2: annotation. No header.')
        parser.add_argument('--noAnnotationTag', metavar='noAnnotationTag', type=str, default = "Other",
                             help='Text to write for the proteins that are not in provided annotation files. Default = "Other"')
#         parser.add_argument('--domainRegex', metavar='domainRegex', type=str, default = "*",
#                              help=' E.g. RRM_*, KH_*, zf-CCCH*, zf-CCHC*, S1, PWI, PUF, SAM_*.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        # init
        run = RBPDomainScore( args.annotationFile, args.catRAPIDFile, args.rainetDB, args.outputFolder, 
                              args.maskMultiple, args.annotationColumn, args.idColumn, args.extraAnnotation, args.noAnnotationTag)
    
        # get cross references from rainet DB
        Timer.get_instance().step( "Read protein cross references..")    
        uniprotACs, protCrossReference = run.protein_cross_references()
    
        # read domain annotations file
        Timer.get_instance().step( "Reading protein annotations..")    
        proteinAnnotation = run.read_domain_annotation_file( protCrossReference)

        # read extra annotations file, if provided
        proteinAnnotation = run.read_extra_annotation_file( proteinAnnotation)      

        # read catrapid file and write output
        Timer.get_instance().step( "Reading catrapid interaction file..")    
        run.read_catrapid_file( proteinAnnotation)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

