from fr.tagc.rainet.core.execution.ExecutionStrategy import ExecutionStrategy
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from sqlalchemy import or_,and_
from fr.tagc.rainet.core.data.DBParameter import DBParameter
from fr.tagc.rainet.core.data.GeneOntology import GeneOntology
from fr.tagc.rainet.core.data.GeneSymbol import GeneSymbol
from fr.tagc.rainet.core.data.KEGGPathway import KEGGPathway
from fr.tagc.rainet.core.data.NetworkModuleAnnotation import NetworkModuleAnnotation
from fr.tagc.rainet.core.data.NetworkModule import NetworkModule
from fr.tagc.rainet.core.data.OCGPartitionAnalysis import OCGPartitionAnalysis
from fr.tagc.rainet.core.data.PartitionAnalysis import PartitionAnalysis
from fr.tagc.rainet.core.data.PPINetworkInteraction import PPINetworkInteraction
from fr.tagc.rainet.core.data.PPINetwork import PPINetwork
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.ProteinDomain import ProteinDomain
from fr.tagc.rainet.core.data.ProteinGOAnnotation import ProteinGOAnnotation
from fr.tagc.rainet.core.data.ProteinInteraction import ProteinInteraction
from fr.tagc.rainet.core.data.ProteinIsoform import ProteinIsoform
from fr.tagc.rainet.core.data.ProteinKEGGAnnotation import ProteinKEGGAnnotation
from fr.tagc.rainet.core.data.ProteinNetworkModule import ProteinNetworkModule
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinReactomeAnnotation import ProteinReactomeAnnotation
from fr.tagc.rainet.core.data.ReactomePathway import ReactomePathway
from fr.tagc.rainet.core.data.SynonymGeneSymbol import SynonymGeneSymbol
from fr.tagc.rainet.core.data.TableStatus import TableStatus




# #
# This class define the Strategy enabling user to perform DB query interactively
class InteractiveQueryStrategy( ExecutionStrategy ):
    

    # #
    # The Strategy execution method
    def execute( self ):
        
        # Get the option values
        self.DBPath = OptionManager.get_instance().get_option( OptionConstants.OPTION_DB_NAME )
        query_file_path = OptionManager.get_instance().get_option( OptionConstants.OPTION_QUERY_FILE )
        species = OptionManager.get_instance().get_option( OptionConstants.OPTION_SPECIES )
        
        # Check if a species has been provided
        if species == None or len( species) == 0:
            raise RainetException( "InteractiveQueryStrategy.execute: You must specify a species in your command line. Please check help.")
        
        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath( self.DBPath)
        sql_session = SQLManager.get_instance().get_session()
        
        # Check if the species in the DB correspond to the species given by the user
        SQLManager.check_species(sql_session, species, True)
        
        # Read the query from file
        Logger.get_instance().info("InteractiveQueryStrategy.execute : Reading query...")
        query_string = self.read_query( query_file_path)
        
        # Get the result of the query on DB
        Logger.get_instance().info("InteractiveQueryStrategy.execute : Querying database...")
        query_result = self.perform_query( sql_session, query_string)
        
        # Export query result to file
        Logger.get_instance().info("InteractiveQueryStrategy.execute : Exporting result...")
        self.export_query_result( query_result, query_file_path)
        
        Logger.get_instance().info("InteractiveQueryStrategy.execute : Finished.")
        
        
    # #
    # Read the query from the file containing its definition
    #
    # @param query_file_path : string - the path to the file with the query definition
    #
    # @return string - The query string if the file was opened and correctly read, None if not
    def read_query(self, query_file_path):
        
        try:
            file_handler = FileUtils.open_text_r( query_file_path)
        except RainetException:
            Logger.get_instance().error( "InteractiveQueryStrategy.read_query : Unable to open query file : " + query_file_path)
            raise RainetException( "InteractiveQueryStrategy.read_query : Unable to open query file : " + query_file_path)
        
        line = file_handler.readline()
        query_string = ''
        while line != None and line != '':
            query_string = query_string + line[0:-1]
            line = file_handler.readline()
            
        return query_string
        
    # #
    # Execute the query on the DB and get the result
    #
    # @param sql_session : SQLSession - the SQL session to the DB
    # @param query_string : string - The query string
    #
    # @return undefined - The result of the query on the DB (should be an array or an integer) 
    #
    # @raise RainetException if the query failed to execute
    def perform_query(self, sql_session, query_string):
        
        full_query = 'sql_session.' + query_string
        Logger.get_instance().info( "InteractiveQueryStrategy.perform_query : query is '" + full_query + "'")
        
        try:
            query_result = eval( full_query)
        except Exception as ex:
            raise RainetException( "InteractiveQueryStrategy.perform_query : Exception occurred during query on DB", ex)
        
        return query_result
    
    # #
    # Export the result in a text file
    #
    # @param query_result : undefined - the query result (may be an array or an integer)
    # @param query_file_path : string - the path to the file containing the query definition
    #
    def export_query_result(self, query_result, query_file_path):
        
        # Check if the query result contains something
        if query_result == None:
            Logger.get_instance().error( "InteractiveQueryStrategy.export_query_result : the query result is None.")
            return
        
        if len( query_result) == 0:
            Logger.get_instance().error( "InteractiveQueryStrategy.export_query_result : the query result is empty.")
            return
        
        # Open the output file
        result_file_path = FileUtils.remove_extension( query_file_path) + "_result.csv"
        try:
            result_file_handler = FileUtils.open_text_w( result_file_path)
        except RainetException:
            Logger.get_instance().error( "InteractiveQueryStrategy.export_query_result : Unable to create result file " + result_file_path)
            return
        
        # Get the list of columns in the queried table
        first_object = query_result[ 0]
        column_list = first_object.__table__.columns
        
        # Write the column names as headers in the output file
        line = ''
        for column in column_list:
            line += column.name + "\t"
        result_file_handler.write( line + "\n")
        
        Logger.get_instance().error( "InteractiveQueryStrategy.export_query_result : query result contains " + str( len( query_result)) + " entries")
        
        # Write the results to the output file
        for obj in query_result:
            line = ''
            for column in column_list:
                line += str( getattr( obj, column.name)) + "\t"
            result_file_handler.write( line + "\n")
        
        #Close the output file
        result_file_handler.close()