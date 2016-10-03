import sys
import os
import argparse


from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil


#===============================================================================
# Started 03-Oct-2016 
# Diogo Ribeiro
DESC_COMMENT = "Wrapper of NetworkScoreAnalysis to run for several top numbers."
SCRIPT_NAME = "networkAnalysisWrapper.py"
#===============================================================================

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
        parser.add_argument('script', metavar='script', type=str,
                             help='Script to run.')
        parser.add_argument('networkFile', metavar='networkFile', type=str,
                             help='PPI Network file to use.')
        parser.add_argument('catrapidFile', metavar='catrapidFile', type=str,
                             help='Catrapid file to use.')
        parser.add_argument('topPartnersList', metavar='topPartnersList', type=str,
                             help='List of top partners parameters to use. Comma-separated E.g. 2,5,10,20 .')
        parser.add_argument('numberRandomizations', metavar='numberRandomizations', type=str,
                             help='Number of randomizations to use.')
           
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        # parse topPartnersList
        topPartnersList = args.topPartnersList.split(",")

        for topPartner in topPartnersList:

            outFolder = "topPartners" + str(topPartner)
            
            # make output folder
            if not os.path.exists( outFolder):
                os.mkdir( outFolder)
            
            cmd = "python %s %s %s %s %s --numberRandomizations %s > %s/log.out 2> %s/log.err &" % ( args.script, \
                    args.networkFile, args.catrapidFile, topPartner, outFolder, args.numberRandomizations, outFolder, outFolder)
            
            print cmd
            SubprocessUtil.run_command( cmd)


        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

