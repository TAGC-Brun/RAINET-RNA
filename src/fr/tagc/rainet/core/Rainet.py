
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.property.PropertyManager import PropertyManager
from fr.tagc.rainet.core.util import Constants

from fr.tagc.rainet.core.execution.InteractiveQueryStrategy import InteractiveQueryStrategy
from fr.tagc.rainet.core.execution.InsertionStrategy import InsertionStrategy


##
# This is the main class of the Rainet project. It contains methods used to insert data to the database
# and make several type of analysis
#
class Rainet( object ):
    

    ##
    # Execute the right strategy according the user command line        
    def execute(self):      
        
        strategy_command = OptionManager.get_instance().get_strategy()
        
        if strategy_command != None:
            try:
                strategy = eval( strategy_command + "Strategy()")
            except Exception:
                raise RainetException( "Rainet.execute : No strategy associated to keyword " + str(strategy_command))    
        else:
            raise RainetException( "Rainet.execute : No strategy was defined: aborting")
    
        try:
            strategy.execute()
        except RainetException as raie:
            Logger.get_instance().error( "Rainet.execute: An exception occurred executing the command:\n" + raie.to_string())
    
#===============================================================================
# The main function
#===============================================================================
if __name__ == '__main__':
        
    try:
        
        # Create Logger instance by using the first log action.
        Logger.get_instance().info( "Rainet : Starting..." )
    
        # Store the options
        OptionManager.get_instance().initialize( )
    
        # Set the level of verbosity
        Logger.get_instance().set_level( OptionManager.get_instance().get_option( OptionConstants.OPTION_VERBOSITY))
    
        # Instantiate Rainet with the correct database path
        rainet = Rainet()
        
        # Insert the data to database
        rainet.execute()
        
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of Rainet. Aborting :\n" + rainet.to_string())
        
    
