from optparse import OptionParser

import OptionConstants

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from os.path import sys
from fr.tagc.rainet.core.util.log.Logger import Logger


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
        self.optionDict = None
        
    # #
    # Initialize the manager with the option parser that contains the option information
    #
    # @param option_parser : OptionParser - the instance of OptionParser built with the command line arguments
    def initialize( self ):
        
        # Get the main keyword that defines the strategy
        self.strategy = sys.argv[1]
        Logger.get_instance().info("Chosen strategy = " + self.strategy)
        
        # If the strategy is not known, raise an exception
        if self.strategy not in OptionConstants.STRATEGIES_LIST:
            raise RainetException( "OptionManager.initialize() : The main keyword is not correct : '" + self.strategy + "'. Should be one of " + str(OptionConstants.STRATEGIES_LIST))
        
        # Build an option parser to collect the option values
        option_parser = OptionParser()
        for current_prop_list in OptionConstants.OPTION_LIST[ self.strategy]:
            option_parser.add_option( current_prop_list[0],
                                          current_prop_list[1],
                                          action = current_prop_list[2],
                                          type = current_prop_list[3],
                                          dest = current_prop_list[4],
                                          default = current_prop_list[5],
                                          help = current_prop_list[6] )
        
        
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
    # @raise RainetException : When the option is not available and not_none is set to True
    def get_option( self, option_name, not_none = False ):
        
        if option_name in self.optionDict.keys():
            option = self.optionDict[ option_name]
            if option != None:
                return option
            
        if not_none == True:
            raise RainetException( "OptionManager.get_option: You must provide the option '" + option_name + "'. See help for more information." )
        
        return None

    # #
    # Set the value of the given parameter in the option dictionary
    #
    # @param option_name : string - The name of the option
    # @param option_value : string - The value of the option
    # 
    def set_option( self, option_name, option_value):
        
        if self.optionDict == None:
            self.optionDict = {}
            
        if option_name != None and len( option_name) > 0:
            if option_value != None:
                self.optionDict[ option_name] = option_value
            else: Logger.get_instance().warning( "OptionManager.set_option: Trying to pass a None or void value on option dictionary for key : " + option_name) 
        else:
            Logger.get_instance().warning( "OptionManager.set_option: Trying to pass a None or void key on option dictionary.")

    # #
    # Returns the singleton instance
    #
    # @return the singleton instance
    @staticmethod
    def get_instance():

        if OptionManager.__instance == None:
            OptionManager.__instance = OptionManager()
        return OptionManager.__instance