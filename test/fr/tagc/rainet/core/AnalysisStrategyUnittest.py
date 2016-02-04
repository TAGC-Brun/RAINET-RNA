import unittest
from fr.tagc.rainet.core.Rainet import Rainet
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.execution.AnalysisStrategy import AnalysisStrategy


class AnalysisStrategyUnittest(unittest.TestCase):
    
    def setUp(self):
        
        # Set the options
        optionManager = OptionManager.get_instance()
        optionManager.set_option(OptionConstants.OPTION_VERBOSITY,"debug")
        optionManager.set_option(OptionConstants.OPTION_DB_NAME,"/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/rainet2016-02-03.human_unittest.sqlite")
        optionManager.set_option(OptionConstants.OPTION_SPECIES,"human")

        # Set the level of verbosity
        Logger.get_instance().set_level( OptionManager.get_instance().get_option( OptionConstants.OPTION_VERBOSITY))

        #create instance of strategy    
        self.strategy = AnalysisStrategy()
        

    def testOne(self):
        
        self.strategy.execute()

    
    def tearDown(self):
        pass
    
    
    