
import unittest
import os
# import pandas as pd
# import glob

# from fr.tagc.rainet.core.execution.processing.catrapid.ReadCatrapid import ReadCatrapid

# #
# Unittesting the ... script on a specific validated dataset. 
#
class TemplateUnittest(unittest.TestCase):

    # Constants with default paramters        
        
    # #
    # Runs before each test
    # name of this function needs forcely to be 'setUp'
    def setUp(self):

        
        #===============================================================================
        # Initialise a normal script
        #===============================================================================

        # Set the options

#         self.catRAPIDFile = "test_input/catRAPID_interactions_test.txt"
#         self.outputFolder = "test_output/"
#         # folder containing expected output files
#         self.expectedFolder = "test_expected/"
        
#        self.run = ReadCatrapid(self.catRAPIDFile, self.outputFolder) 

        #===============================================================================
        # Initialise a Rainet.py script
        #===============================================================================

#         optionManager = OptionManager.get_instance()
#         optionManager.set_option(OptionConstants.OPTION_VERBOSITY, "debug")
#         optionManager.set_option(OptionConstants.OPTION_DB_NAME, "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/rainet_testing_DB.sqlite") #rainet2016-06-17.human_expression_wPRI.sqlite")
#         # Set the level of verbosity
#         Logger.get_instance().set_level(OptionManager.get_instance().get_option(OptionConstants.OPTION_VERBOSITY))
# 
#         # Setting up SQL manager
#         SQLManager.get_instance().set_DBpath(OptionManager.get_instance().get_option(OptionConstants.OPTION_DB_NAME))
#         self.sql_session = SQLManager.get_instance().get_session()
#        # create instance of strategy    
#         self.run = EnrichmentAnalysisStrategy()
#         self.run.sql_session = SQLManager.get_instance().get_session()
 
        pass

    # #
    def test_default_params(self):

        print "| test_default_params | "

#         self.run.run( )
                
#         self.assertTrue( round(proteinInteractionsMean["Q7Z419"],1) == 18.0, "assert that..")


    # #
    def test_params_one(self):

        print "| test_params_one | "

#         # check if output files are correct
#         with open(self.outputFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME, "r") as out:                
#             with open(self.expectedFolder + ReadCatrapid.STORED_INTERACTIONS_FILENAME, "r") as exp:
#                 self.assertTrue(out.read() == exp.read(), "assert if report file is correct, by expected content comparison" )
# 
#         with open(self.outputFolder + ReadCatrapid.PROTEIN_INTERACTIONS_FILENAME, "r") as out:                
#             with open(self.expectedFolder + ReadCatrapid.PROTEIN_INTERACTIONS_FILENAME, "r") as exp:
#                 self.assertTrue(out.read() == exp.read(), "assert if report file is correct, by expected content comparison" )
# 
#         # cp /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_output/* /home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/processing/catrapid/test_expected/

        pass

    # #
    def test_report_two(self):

        print "| test_report_two | "

#         optionManager = OptionManager.get_instance()        
#         optionManager.set_option(OptionConstants.OPTION_MINIMUM_INTERACTION_SCORE, "28")
# 
#         self.run.run()
#         
#         # RNA numbers report
#         table = pd.read_table( self.outputFolder + AnalysisStrategy.REPORT_RNA_NUMBERS)
#                 
#         self.assertTrue( table["Total_Genes"][0] == 198, "assert if number of Genes before filter is correct")
#         self.assertTrue( table["Total_Genes"][1] == 8, "assert if number of Genes after filter is correct")        

    def test_hypergeometric_test(self):
        
        print "| test_hypergeometric_test | "

#         x, m, n, k = 5, 10, 50, 10
# 
#         res = "%.5f" % self.run.hypergeometric_test( x, m, n, k)
# 
#         self.assertTrue( res == "0.00067")

    
    # #
    # Runs after each test
    def tearDown(self):
                       
        # Wipe output folder
        cmd = "rm %s/*" % self.outputFolder
        os.system(cmd)
            
       


