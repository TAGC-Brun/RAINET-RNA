import unittest
from fr.tagc.rainet.core.Rainet import Rainet
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.execution.AnalysisStrategy import AnalysisStrategy


class RainetUnittest(unittest.TestCase):
    
    def setUp(self):
        
        # Create Logger instance by using the first log action.
        Logger.get_instance().info( "Rainet : Starting..." )
    
        # Store the options
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_VERBOSITY,"debug")

        #create instance of strategy
    
        strategy = AnalysisStrategy()
#        strategy.execute()

        # Set the level of verbosity
        Logger.get_instance().set_level( OptionManager.get_instance().get_option( OptionConstants.OPTION_VERBOSITY))

        print (OptionManager.get_instance().get_option( OptionConstants.OPTION_LIST))

    def testOne(self):
        
        pass
    
    def tearDown(self):
        pass
    
    
    