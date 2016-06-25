import sys
import os
import argparse

# import numpy as np
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
DESC_COMMENT = "Script to ..."
SCRIPT_NAME = "name_of_script.py"
#===============================================================================

#===============================================================================
# General plan:
# 1)
# 2)
#===============================================================================

#===============================================================================
# Processing notes:
# 1)
# 2)
#===============================================================================

class TemplateScriptWithClass(object):
    
    #===============================================================================
    # TemplateScriptWIthClass Constants
    #===============================================================================

    # significance value 
    SIGN_VALUE = 0.05
        
    #===================================================================
    # Data Manager object Keywords
    #===================================================================

    # Protein / RNA with interaction data
    PRI_PROT_KW = "interactingProteins" # Stores all Proteins in interactions
    PRI_RNA_KW = "interactingRNAs" # Stores RNAs in interactions

    #===================================================================
    # Report files constants       
    #===================================================================

    PARAMETERS_LOG = "parameters.log"

    # Annotation report
    REPORT_PROT_PER_ANNOTATION = "prot_per_annotation.tsv"

       
    def __init__(self, arg_one, arg_two):

        self.argOne = arg_one
        self.argTwo = arg_two

#         # Build a SQL session to DB
#         SQLManager.get_instance().set_DBpath(self.rainetDB)
#         self.sql_session = SQLManager.get_instance().get_session()

#         # make output folder
#         if not os.path.exists( self.outputFolder):
#             os.mkdir( self.outputFolder)


    # #
    # Write description of function
    def function_one( self):
        
        
#         # Query all UniProtAC in database
#         query = self.sql_session.query( Protein.uniprotAC ).all()        

#         with open( self.argOne, "r") as inFile:

        if 1:
            raise RainetException("function_one: error message: ",  TemplateScriptWithClass.PARAMETERS_LOG)


        return TemplateScriptWithClass.PARAMETERS_LOG
   
    # #
    # Central function to run functions in order
    def run( self):

        Logger.get_instance().info( "TemplateScriptWithClass.run: Starting..." )

        #===================================================================
        # Initialising datasets
        #===================================================================

        Timer.get_instance().step( "Initialising interaction datasets.." )        
          
        self.function_one()
         
#         # Container for results of hypergeom tests so that they don't have to be repeated. 
#         Logger.get_instance().info( "EnrichmentAnalysisStrategy.analysis : Number of unique tests performed: %s" % ( len( self.testContainer)) )
   

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
        parser.add_argument('sysArgFile', metavar='sysArgFile', type=str,
                             help='Description.')
        # optional args
        parser.add_argument('--sysArgOptional', metavar='sysArgOptional', type=int, default = 1,
                             help='Description (Default = 1).')
           
        # gets the arguments
        args = parser.parse_args( ) 
    
        # Initialise class
        templateScriptWIthClass = TemplateScriptWithClass( args.sysArgFile, args.sysArgOptional)
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        # run function to..
        Timer.get_instance().step( "Read bogus file..")            
        templateScriptWIthClass.run( )

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + rainet.to_string())

