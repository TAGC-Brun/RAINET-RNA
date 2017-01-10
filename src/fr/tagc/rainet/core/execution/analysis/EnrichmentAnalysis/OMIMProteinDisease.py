
import sys
import os
import argparse

# import numpy as np
# import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
import cmd
# from fr.tagc.rainet.core.util.data.DataManager import DataManager


#===============================================================================
# Started 05-January-2017
# Diogo Ribeiro
DESC_COMMENT = "Script to retrive OMIM data from biomart, get corresponding descriptions and create protein-disease file."
SCRIPT_NAME = "OMIMProteinDisease.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) Run R script to access Biomart and retrieve OMIM-protein correspondence
# 2) Read OMIM mimTitles.txt file to connect OMIM disease ID to disease description
# 3) Return file with protein-disease correspondence
#===============================================================================

#===============================================================================
# Processing notes:
# 1) Since some uniprotAC-disease correspondences are missing, we also use HGNC-disease correspondences
#===============================================================================

class OMIMProteinDisease(object):

    # Constants

    R_SCRIPT_OUTPUT_FILE_NAME = "/OMIM_biomart.txt"
    
    OUTPUT_FILE = "/protein_disease_description.txt"
        
    def __init__(self, mimTitlesFile, outputFolder, rscriptPath):

        self.mimTitlesFile = mimTitlesFile
        self.outputFolder = outputFolder
        self.rscriptPath = rscriptPath

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)

        self.rscriptOutputFile = self.outputFolder + OMIMProteinDisease.R_SCRIPT_OUTPUT_FILE_NAME

        # remove R output file in case it exists, so that it does not overwrite
        if os.path.exists( self.rscriptOutputFile):
            cmd = "rm %s" % self.rscriptOutputFile
            SubprocessUtil.run_command( cmd, verbose = 0)
            
            
    # #
    # Run R script to retrieve OMIM biomart data into a predefined file
    def run_biomart_rscript(self):
        
        cmd = "Rscript %s %s" % ( self.rscriptPath, self.rscriptOutputFile)
        
        returnCode = SubprocessUtil.run_command( cmd, verbose = 1)
        
        print "run_biomart_rscript: return code %s" % returnCode

    # #
    # Read output file from R script biomart
    def read_biomart_data(self):

        #===============================================================================
        # R script biomart file
        #===============================================================================
        # Example format
        # ensembl_gene_id mim_morbid      hgnc_symbol     uniprot_swissprot
        # ENSG00000281766 NA      RYBP    
        # ENSG00000281518 NA      FOXO6  
        # ENSG00000281917 616095  SLC16A1 P53985

        # Note: there is one entry for each Ensembl ID, these may have corresponding OMIM and UniprotAC or not.
        # We are interested in getting only correspondence from uniprotAC to mim_morbid
        # However, there are correspondences between HGNC IDs and mim_morbid, that are not present as uniprotAC, even if the HGNC and uniprotAC have correspondece.
        # Therefore we read both protein IDs, and map them to diseases in a greedy way.

        #===============================================================================
        # First read to get HGNC-uniprotAC correspondence
        #===============================================================================

        hgnc2Uniprot = {} # key -> HGNC ID, value -> uniprotAC
        uniprotACs = set()
        # Note: a single hgncID can have several uniprotACs mapped to it.
        with open( self.rscriptOutputFile, "r") as inFile:
            
            header = inFile.readline()
            
            for line in inFile:
                spl = line.strip().split( "\t")
                
                if len( spl) < 4:
                    # incomplete line, no uniprotAC, or no HGNC and uniprot
                    continue
                
                hgncID = spl[2]
                uniprotAC = spl[3]

                if len( hgncID) == 0:
                    # if there is uniprot entry but no HGNC entry
                    continue

                if hgncID not in hgnc2Uniprot:
                    hgnc2Uniprot[ hgncID] = set()

                hgnc2Uniprot[ hgncID].add( uniprotAC)
                
                uniprotACs.add( uniprotAC)

        print "read_biomart_data: %s hgnc IDs read (with uniprotAC correspondence)" % len( hgnc2Uniprot)
        print "read_biomart_data: %s uniprotACs read" % len( uniprotACs)

        #===============================================================================
        # Second read to get disease-uniprotAC correspondence
        #===============================================================================
        # Note: using HGNC ID correspondence between uniprotAC and HGNC ID

        proteinOMIMDict = {} # key -> uniprotAC, value -> OMIM number / ID
        omimIDs = set()
        nLines = 0
        
        with open( self.rscriptOutputFile, "r") as inFile:
            
            header = inFile.readline()
            nLines += 1
            
            for line in inFile:
                spl = line.strip().split( "\t")

                uniprotAC = ""
                hgncID = ""
                omimID = ""

                nLines += 1
                
                if len( spl) < 3:
                    # incomplete line, no uniprotAC and no hgnc
                    continue
                
                omimID = spl[1]

                if omimID == "NA":
                    # nothing to retrieve
                    continue        

                if len( spl) > 3:
                    # if there is uniprotAC, go for direct correspondence
                    uniprotAC = spl[3]
  
                    if uniprotAC not in proteinOMIMDict:
                        proteinOMIMDict[ uniprotAC] = set()
     
                    proteinOMIMDict[ uniprotAC].add( omimID)
                    omimIDs.add( omimID)

                ## Greedy approach
                # regardless of uniprotAC, make HGNC-disease correspondence, and add it as uniprotAC-disease by hgnc-uniprotac correspondence
                # a posteriori we realised that this approach only brought one extra uniprotac-disease association

                hgncID = spl[2]
 
                if len( hgncID) > 0:
                    # if there is hgnc to uniprotac correspondence, use that entry as well
                    if hgncID in hgnc2Uniprot:
                        uniprots = hgnc2Uniprot[ hgncID]
                        for unip in uniprots:
                            if unip not in proteinOMIMDict:
                                proteinOMIMDict[ unip] = set()
                            proteinOMIMDict[ unip].add( omimID)
                            omimIDs.add( omimID)
                    else:
                        # this should be only proteins that dont have any reviewed uniprot entry (but may have unreviewed)
                        pass

        
        print "read_biomart_data: %s uniprotACs in OMIM" % len( proteinOMIMDict)
        print "read_biomart_data: %s OMIM IDs" % len( omimIDs)

        self.proteinOMIMDict = proteinOMIMDict


    # #
    # Read mimTitles.txt file from OMIM which connects OMIM ID to a disease (and other) description
    def read_mim_titles_file(self):

        #===============================================================================
        # Read mimTitles.txt
        #===============================================================================
        # Example file        
        # # Copyright (c) 1966-2016 Johns Hopkins University. Use of this file adheres to the terms specified at http://omim.org/help/agreement.
        # # Generated: 2016-06-29
        # # Prefix        Mim Number      Preferred Title; symbol Alternative Title(s); symbol(s) Included Title(s); symbols
        # NULL    100050  AARSKOG SYNDROME, AUTOSOMAL DOMINANT            
        # Percent 100070  AORTIC ANEURYSM, FAMILIAL ABDOMINAL, 1; AAA1    ANEURYSM, ABDOMINAL AORTIC; AAA;; ABDOMINAL AORTIC ANEURYSM     
        # Number Sign     100100  PRUNE BELLY SYNDROME; PBS       ABDOMINAL MUSCLES, ABSENCE OF, WITH URINARY TRACT ABNORMALITY AND CRYPTORCHIDISM;; EAGLE-BARRETT SYNDROME; EGBRS        
        # NULL    100200  ABDUCENS PALSY          

        # Note: comment lines start with "#" and occur in begginning and end of file

        minTitlesDict = {} # key -> OMIM ID, value -> description
        
        with open( self.mimTitlesFile, "r") as inFile:
            
            for line in inFile:

                # reset descriptions
                desc = ""
                altDesc = ""
                incDesc = ""

                # skip comment lines
                if line.startswith("#"):
                    continue

                spl = line.strip().split("\t")
                
                mimNumber = spl[1]
                desc = spl[2]
                if len(spl) > 3:
                    altDesc = spl[3]
                if len(spl) > 4:
                    incDesc = spl[4]

                # merge original description, alternative and included
                descText = desc + " | " + altDesc + " | " + incDesc

                minTitlesDict[ mimNumber] = descText

        print "read_mim_titles_file: %s OMIM ID entries read." % len( minTitlesDict)

        self.minTitlesDict = minTitlesDict


    # #
    # Write output file with uniprotACs and their associated diseases
    def write_output(self):
        
        #===============================================================================
        # Output file
        #===============================================================================
        # proteinID\tOMIM_ID\tOMIM_description        
        
        # Note: output file only contains proteins that have any associated disease
        # Note: a protein may have several diseases, the output file contains one line per disease of a protein
        
        outFile = open( self.outputFolder + OMIMProteinDisease.OUTPUT_FILE, "w")
        
        notFound = set()

        nEntries = 0
        
        for prot in self.proteinOMIMDict:
            for omimID in self.proteinOMIMDict[ prot]:
                
                # Important: all OMIM IDs should be found if using updated mimTitles.txt file
                if omimID in self.minTitlesDict:
                    outFile.write( "%s\t%s\t%s\n" % (prot, omimID, self.minTitlesDict[ omimID]) )
                    nEntries += 1
                else:
                    Logger.get_instance().warning( "write_output: %s OMIM ID not found in mimTitles.txt" % omimID)

                    notFound.add( omimID)

        print "write_output: %s OMIM IDs not found." % len( notFound)
        print "write_output: wrote %s protein-disease association entries." % nEntries

        outFile.close()


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
        parser.add_argument('mimTitlesFile', metavar='mimTitlesFile', type=str,
                             help='Path to mimTitles.txt file downloadable from OMIM database.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Output folder.')
        parser.add_argument('--rscriptPath', metavar='rscriptPath', type=str, default = "OMIM_biomart.R",
                             help='Path to R script from retrieving biomart OMIM data.')           
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        instance = OMIMProteinDisease( args.mimTitlesFile, args.outputFolder, args.rscriptPath)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        Timer.get_instance().step( "Retrieve biomart data..")
        instance.run_biomart_rscript( )

        Timer.get_instance().step( "Read biomart data..")
        instance.read_biomart_data( )

        Timer.get_instance().step( "Read mimTitles OMIM file..")            
        instance.read_mim_titles_file()

        Timer.get_instance().step( "Write output file..")            
        instance.write_output()
        
        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

