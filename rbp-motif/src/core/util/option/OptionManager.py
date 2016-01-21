from optparse import OptionParser

import OptionConstants

from core.util.exception.RbpmotifException import RbpmotifException
from os.path import sys
from core.util.log.Logger import Logger


# # This class is a singleton aiming to manage the execution options provided by the user command.
# This class is in fact a wrapper for the built-in OptionParser used to parse the command line
#
class OptionManager( object ) :

    __instance = None
    
    # #
    #
    def __init__( self ):
        
        self.optionParser = None
        self.args = None
        
    # #
    # Initialize the manager with the option parser that contains the option information
    #
    # @param option_parser : OptionParser - the instance of OptionParser built with the command line arguments
    def initialize( self ):
        
        # Build an option parser to collect the option values
        option_parser = OptionParser()
        for current_prop_list in OptionConstants.OPTION_LIST:
            option_parser.add_option( current_prop_list[0],
                                          current_prop_list[1],
                                          action = current_prop_list[2],
                                          type = current_prop_list[3],
                                          dest = current_prop_list[4],
                                          default = current_prop_list[5],
                                          help = current_prop_list[6],
                                          nargs = current_prop_list[7])
        
        
        # Get the various option values into a dictionary
        ( opts, args ) = option_parser.parse_args()
        self.optionDict = vars( opts )
        
        
    ## 
    # Returns the value of 
    #
    def get_strategy(self):
        
        return self.strategy
    
    # #
    # Returns the value of the option with the provided option_name.
    # if the option is not available, return None or an Exception
    #
    # @param option_name : string - The name of the option to get
    # @param not_none : boolean - (optional) indicates if a None value can be returned when
    #                     option is not available. If no, an exception is raised
    # @return The value of the option if present.
    # @raise RbpmotifException : When the option is not available and not_none is set to True
    def get_option( self, option_name, not_none = False ):
        
        if option_name in self.optionDict.keys():
            option = self.optionDict[ option_name]
            if option != None:
                return option
            
        if not_none == True:
            raise RbpmotifException( "OptionManager.get_option: You must provide the option '" + option_name + "'. See help for more information." )
        
        return None

    # #
    # Returns the singleton instance
    #
    # @return the singleton instance
    @staticmethod
    def get_instance():

        if OptionManager.__instance == None:
            OptionManager.__instance = OptionManager()
        return OptionManager.__instance
