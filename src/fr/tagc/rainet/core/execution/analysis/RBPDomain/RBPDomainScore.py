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
DESC_COMMENT = "Script to attribute nucleic-acid binding domains to proteins of a processed catrapid file."
SCRIPT_NAME = "RBPDomainScore.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Read file with protein domain annotation
# 2) Read processed catRAPID file, add annotation
#===============================================================================

class RBPDomainScore(object):
    
    # database for which we want cross references
    XREF_SOURCE_DB = "Ensembl_PRO"
    UNIPROT_OUTPUT_FILE = "/proteins_uniprotac.tsv"
    ANNOTATION_OUTPUT_FILE = "/annotated_interactions.tsv"
    NO_ANNOTATION_TAG = "Non-RBP"
    SEVERAL_ANNOTATION_TAG = "several_annotations"
       
    def __init__(self, domain_annotation_file, catrapid_file, rainet_db, output_folder, mask_multiple):

        self.domainAnnotationFile = domain_annotation_file
        self.catRAPIDFile = catrapid_file
        self.rainetDB = rainet_db
        self.outputFolder = output_folder
        self.maskMultiple = mask_multiple

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()


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

        with open( self.domainAnnotationFile, "r") as inFile:
            # skip header
            inFile.readline()

            notFound = set()
            
            for line in inFile:
                line = line.strip()
                lineCounter+=1
                
                geneID,pepID,pfamIDs,pfamNames,domainClass = line.split( "\t")
                                
                if not pepID.startswith("ENSP"):
                    raise RainetException("read_domain_annotation_file: Peptide ID is incorrect:", pepID)

                # retrieve corresponding uniprotACs
                if pepID in protCrossReference:
                    protIDs = protCrossReference[ pepID]
                else:
                    notFound.add( pepID)
                    continue

                # there can be several domain classes for the same protein
                domainClasses = domainClass.split( ",")

                outFile.write( "%s\t%s\t%s\n" % ( pepID, ",".join( protIDs), domainClass ))

                for cla in domainClasses:
                    # there can be several protein uniprotAC for a single Ensembl peptideID
                    for protID in protIDs:
                        if protID not in proteinAnnotation:
                            proteinAnnotation[ protID] = set()
                        proteinAnnotation[ protID].add( cla)

            print "number of entries read:", lineCounter
            print "number of peptide IDs not found: ", len(notFound)
            print "number of proteins with annotation:", len(proteinAnnotation)

        outFile.close()

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
        # output file header
        outFile.write( "uniprotac\tscore\tannotation\n")
          
        #=======================================================================
        # read file
        #=======================================================================

        with open( self.catRAPIDFile, "r") as inFile:
            # skip header
            inFile.readline()
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
                        outputText += "%s\t%s\t%s\n" % (protID, score, annotation)

                    # if wanting information for all annotations
                    elif len( proteinAnnotation[ protID]) > 0:
                        # get all annotations
                        annotations = ",".join( list(proteinAnnotation[ protID]) )
                        outputText =  "%s\t%s\t%s\n" % (protID, score, annotations)
                    else:
                        raise RainetException("read_catrapid_file: Annotation information is incorrect. ", proteinAnnotation[ protID])
                # if there is no annotation for protein
                else:
                    annotation = RBPDomainScore.NO_ANNOTATION_TAG
                    outputText+= "%s\t%s\t%s\n" % (protID, score, annotation)
                
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
        parser.add_argument('domainAnnotationFile', metavar='domainAnnotationFile', type=str,
                             help='TSV file with domain annotation per protein. E.g. gene id prot id pfam id domain name     domain class.')
        parser.add_argument('catRAPIDFile', metavar='catRAPIDFile', type=str,
                             help='proteinInteractions.tsv output file from ReadCatrapid.py. E.g. uniprotac\tmean_score')
        parser.add_argument('rainetDB', metavar='rainetDB', type=str, help='Path to RAINET database to be used.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
        parser.add_argument('--maskMultiple', metavar='maskMultiple', type=int, default = 1,
                             help='Whether to mask protein annotations when having more than one annotation (val = 1), or display all annotations separated by comma (val = 0). (default = 1).')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        # init
        run = RBPDomainScore( args.domainAnnotationFile, args.catRAPIDFile, args.rainetDB, args.outputFolder, args.maskMultiple)
    
        # get cross references from rainet DB
        Timer.get_instance().step( "Read protein cross references..")    
        uniprotACs, protCrossReference = run.protein_cross_references()
    
        # read domain annotations file
        Timer.get_instance().step( "Reading protein domain annotations..")    
        proteinAnnotation = run.read_domain_annotation_file( protCrossReference)

        # read catrapid file and write output
        Timer.get_instance().step( "reading catrapid interaction file..")    
        run.read_catrapid_file( proteinAnnotation)

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

