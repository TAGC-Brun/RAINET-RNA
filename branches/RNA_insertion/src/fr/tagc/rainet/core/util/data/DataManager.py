from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data import DataConstants

from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.RNACrossReference import RNACrossReference
from fr.tagc.rainet.core.data.MRNA import MRNA
from fr.tagc.rainet.core.data.LncRNA import LncRNA
from fr.tagc.rainet.core.data.OtherRNA import OtherRNA
from fr.tagc.rainet.core.data.ProteinRNAInteractionCatRAPID import ProteinRNAInteractionCatRAPID
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference


# #
# This class is a singleton aiming to manage the data retrieve from internal database queries.
# e.g. storing in memory a list of RNAs for easier accession.
class DataManager( object ) :

    __instance = None
    
    # #
    # Constructor of the Factory
    #
    def __init__( self ):

        #keys of dictionary, given by user, will point to query result
        self.data = {} 
        

    # #
    # Make SQL query to database and store query results.
    #
    # @param query_string: string - The query string
    #
    # @raise RainetException if the query failed to execute
    def perform_query(self, keyword, query_string):

        sql_session = SQLManager.get_instance().get_session()

        full_query = 'sql_session.' + query_string
        
        Logger.get_instance().info( "DataManager.init : query is '" + full_query + "'")

        try:
            query_result = eval( full_query)
        except Exception as ex:
            raise RainetException( "DataManager.init : Exception occurred during query on DB", ex)

        self.data[keyword] = query_result

    # #
    # Process previously stored results into a key-value dictionary of lists
    #
    # @param keyword : String - the data dictionary keyword to access the data
    # @param col1 : Int - 0-based column index to be used as key in the dict
    # @param col2 : Int - 0-based column index to be used as value in the dict
    #
    # @raise RainetException if the keyword is not present
    def query_to_dict(self, keyword, col1, col2):

        outDict = {}
        
        if keyword in self.data:
            data = self.data[keyword]
        else:
            raise RainetException("DataManager.query_to_dict : Data keyword does not exist"+keyword)            


        for entry in data:
            if col1 > len(entry) or col2 > len(entry):
                raise RainetException("DataManager.query_to_dict : requested columns out of boundaries of query result: "+entry)            
            
            keyItem = str(entry[col1])
            valueItem = str(entry[col2])

            if keyItem not in outDict:
                outDict[keyItem] = []
            outDict[keyItem].append(valueItem)

        self.data[keyword] = outDict

    # #
    # Get previously stored results
    #
    # @param keyword: the data dictionary keyword to access the data
    #
    # @raise RainetException if the keyword is not present
    def get_data(self, keyword):
        
        try:
            keyword = str(keyword)
        except:
            raise RainetException("DataManager.get_data : keyword must be a string.")
        
        if keyword not in self.data:
            raise RainetException("DataManager.get_data : Data keyword does not exist"+keyword)
            
        return self.data[keyword]

    # #
    # Delete previously stored results
    #
    # @param keyword: the data dictionary keyword to access the data
    #
    # @raise RainetException if the keyword is not present
    def delete_data(self, keyword):
        
        try:
            keyword = str(keyword)
        except:
            raise RainetException("DataManager.delete_data : keyword must be a string.")
        
        if keyword in self.data:
            del self.data[keyword]
        else:
            raise RainetException("DataManager.delete_data : Data keyword does not exist"+keyword)
            

    # #
    # Returns the singleton instance
    #
    # @return the singleton instance        
    @staticmethod
    def get_instance():

        if DataManager.__instance == None:
            DataManager.__instance = DataManager()
        return DataManager.__instance

