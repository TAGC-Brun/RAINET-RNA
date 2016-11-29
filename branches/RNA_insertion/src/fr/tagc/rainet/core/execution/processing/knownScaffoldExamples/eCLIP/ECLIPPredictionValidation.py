import argparse
import pandas as pd

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

#===============================================================================
# Started 29-Nov-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to see if catRAPID predictions distinguish eCLIP interactions."
SCRIPT_NAME = "ECLIPPredictionValidation.py"
#===============================================================================

#===============================================================================
# General plan:
# 1) 
#===============================================================================
#===============================================================================
# Processing notes:
# 
#===============================================================================


class ECLIPPredictionValidation( object ):
    
    def __init__(self, catrapidFile, eClipFile, outputFolder):

        self.catrapidFile = catrapidFile
        self.eClipFile = eClipFile
        self.outputFolder = outputFolder

        
    # #
    # Read eClip file and retrieve list of proteins interacting with RNAs
    def read_eclip_file(self):
           
        # Example input
        # ENST00000169298|O00425  10

        interactingPairs = {} # key -> pair of transcriptID and proteinID, val -> count of interactions
        proteinSet = set()
        transcriptSet = set()

        with open( self.eClipFile, "r") as inFile:
            for line in inFile:
                spl = line.strip().split( "\t")
                spl2 = spl[0].split( "|")
                transcriptID = spl2[ 0]
                proteinID = spl2[ 1]
                
                proteinSet.add(proteinID)
                transcriptSet.add( transcriptID)
                
                pair = transcriptID + "|" + proteinID
                
                interactingPairs[ pair] = spl[1]
  
        print "read_eclip_file: Total number of interacting pairs:",len(interactingPairs)
        print "read_eclip_file: Total number of interacting RNAs:",len(transcriptSet)
        print "read_eclip_file: Total number of interacting proteins:",len(proteinSet)
        
        self.eclipPairs = interactingPairs


    # #
    # Read catRAPID file.
    # Updated for new catRAPID format. No need for cross references.
    def read_catrapid_file(self):

        # output file
        outFile = open( run.outputFolder + "/scores.tsv", "w")
        outFile.write( "pairID\tcatrapid_score\tin_validated_set\n")

        # E.g.: sp|Q6P6C2|ALKB5_HUMAN ENST00000559683   47.85   0.93    0.23

        proteinSet = set()
        transcriptSet = set()

        countLines = 0
        
        countInValidated = 0
        
        with open( self.catrapidFile, "r") as f:
            for line in f:
                spl = line.split(" ")

                countLines+= 1 

                if countLines % 10000000 == 0:
                    print "Processed %s interactions" % countLines
                               
                proteinID = spl[0].split( "|")[1]
                spl2 = spl[1].split( "\t")
                transcriptID = spl2[0]
                intScore = float( spl2[1])
                
                pair = transcriptID + "|" + proteinID

                proteinSet.add(proteinID)
                transcriptSet.add( transcriptID)

                if pair in self.eclipPairs:
                    inValidated = 1
                    countInValidated += 1
                else:
                    inValidated = 0
   
                outFile.write( "%s\t%s\t%s\n" % ( pair, intScore, inValidated) )
             
        outFile.close()

        print "read_catrapid_file: Number of protein-RNA pairs in catRAPID: ", len( countLines)
        print "read_catrapid_file: Number of proteins: ", len( proteinSet)
        print "read_catrapid_file: Number of transcripts: ", len( transcriptSet)
        print "read_catrapid_file: Number of interactions with experimental data: ", countInValidated



if __name__ == "__main__":
    
    try:
        # Create Logger instance by using the first log action.
        Logger.get_instance().info( " Starting..." )

        #===============================================================================
        # Get input arguments, initialise class
        #===============================================================================
        parser = argparse.ArgumentParser(description='# ') 

        # positional args
        parser.add_argument('catRAPIDFile', metavar='catRAPIDFile', type=str,
                             help='File path of CatRAPID omics/fragments results from the webserver.')
        parser.add_argument('eClipFile', metavar='eClipFile', type=str,
                             help='File path of eClip file from "process_all_eclip_files.py".')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Folder where to write output files.')
        
        #gets the arguments
        args = parser.parse_args( ) 

        # Initialise class
        run = ECLIPPredictionValidation( args.catRAPIDFile, args.eClipFile, args.outputFolder)

        #===============================================================================
        # Run analysis / processing
        #===============================================================================
         
        # Start chrono
        Timer.get_instance().start_chrono()
 
        Timer.get_instance().step( "reading eCLIP file..")    

        run.read_eclip_file()

        Timer.get_instance().step( "reading catRAPID file..")    

        run.read_catrapid_file()


    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution. Aborting :\n" + rainet.to_string())

    # Stop the chrono      
    Timer.get_instance().stop_chrono( " Finished" )

