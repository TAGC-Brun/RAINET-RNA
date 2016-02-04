
from fr.tagc.rainet.core.util.file.FileUtils import FileUtils 
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.factory.DataFactory import DataFactory
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util import Constants
from fr.tagc.rainet.core.util.time.Timer import Timer

## This class is a parser of tab separated text file.
# The class contains methods allowing to parse any TSV file by creating a specific type of Object
# for each line found
# The parser have to received as input:
# - The path to the file to parse
# - The information on the presence or not of column headers in the file
# - The ordered list of column headers (if headers are not present or if user want to overwrite file header values)
# - The name of the Class corresponding to Object to be created for each file line
# - The ordered list of header names to use to call the Class constructor
# - The ordered list of optional values to use to call the Class constructor
# - The comment symbol
# - clean_table : boolean  - The information indicating if the related db table must be cleaned before insertion or not
class TSVParser(object):

    ##
    # This method parses the provided TSV file, looking for headers if required and creating one Object of the requested Class per line
    #
    # @param file_path : string - The path to the TSV file to parse
    # @param has_header : boolean - Indicates if the file to parse has a header line
    # @param header_list : list<string> - Ordered list of headers name (could be None)
    # @param class_name : string - name of the class to instantiate
    # @param parameter_name_list : list<string> - Ordered list of parameters used to instantiate the Class
    # @param parameter_value_list : list<string> - Ordered list of optional parameter values used to instantiate the Class
    # @param comment_symbol : string - Symbol used to comment lines in the file to parse
    # @param clean_table : boolean  - The information indicating if the related db table must be cleaned before insertion or not
    #
    # @return None
    @staticmethod
    def parse_file( file_path, has_headers, group_list, class_name, parameter_name_list, parameter_value_list, comment_symbol, clean_table = False):
        
        # Initialize the status of the insertion
        status = Constants.STATUS_OK
        
        # Open the file to parse in 'read' mode
        input_file = FileUtils.open_text_r( file_path)
        
        # If required, look for the header line (first non-empty line)
        if has_headers:
            line = ''
            while line == '':
                line = input_file.readline()
            # If no headers are provided by the user, use the file headers
            if group_list == None:
                group_list = line.split( "\t")
        
        # Build the map of parameter name to header index
        # If a parameter name is not in the header the index is set to -1
        # and a warning message is sent
        parameter_to_header_index_map = {}
        for parameter_name in parameter_name_list:
            try:
                index = group_list.index( parameter_name)
            except ValueError:
                index = -1
                Logger.get_instance().warning( "TSVParser.parse_file : The parameter '" + parameter_name + "' is not in the header list : " + str(group_list) + ".\n Check whether it is normal or not.")
                status = Constants.STATUS_WARNING
            parameter_to_header_index_map[ parameter_name] = index
                        
        
        # Build the Factory to be used to create the objects
        data_factory = DataFactory( class_name)
        if clean_table:
            data_factory.clean_table()
        
        # Parse the lines in the files and create one object of the requested type per line
        instance_list = []
        line = input_file.readline()
        counter = 0
        while line != None and line != '':
            # Ignore the empty and comment lines
            if line != '' and not line.startswith( comment_symbol):
                line_value_list = line.split( "\t")
                #Check if the line has as many columns as headers
                if len( line_value_list) > len( group_list):
                    raise RainetException( "TSVParser.parse_file : Bad line format in file. Line has not the right number of values.\n"
                        +" Number of values = " + str( len( line_value_list)) + " != Number of headers =" + str( len( group_list)) + ".\n"
                        +" Line =" + line)
                
                #Build the ordered list of value required for Object creation
                #If no parameter_value_list is provided, create a new 'None' one
                # If a parameter_value_list is provided, override only the values
                # corresponding to a parameter found in the file header
                if( parameter_value_list == None):
                    parameter_value_list = [None]*len( parameter_name_list)
                value_index = 0
                for parameter_name in parameter_name_list:
                    parameter_index = parameter_to_header_index_map[ parameter_name]
                    if parameter_index >=0 :
                        parameter_value_list[ value_index] = line_value_list[ parameter_index].strip()
                    value_index = value_index + 1
                    
                #Call for the Object creation
                new_instance = data_factory.create_object_from_tsv( parameter_value_list)
                if new_instance != None:
                    instance_list.append( new_instance)
            
            # Read a new line
            line = input_file.readline()

            
            #dr testing
            counter += 1
            if counter % 100000 == 0:
                Logger.get_instance().info( "TSVParser.parse_file : %s Read %i lines.." % (file_path,counter) )
                Timer.get_instance().time_point()
                SQLManager.get_instance().commit()
                instance_list = []


        # Close the input file
        input_file.close()
        
        # Commit the SQLAlchemy session
        Logger.get_instance().info( "TSVParser.parse_file : Committing SQL session.")
        SQLManager.get_instance().commit();
        
        # Return the list of created instances
        return status
        
        
