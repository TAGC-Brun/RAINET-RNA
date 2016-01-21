
import re

from fr.tagc.rainet.core.util.file.FileUtils import FileUtils 
from fr.tagc.rainet.core.util.factory.DataFactory import DataFactory
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util import Constants

## This class is a parser of tab separated text file.
# The class contains methods allowing to parse any TSV file by creating a specific type of Object
# for each line found
# The parser have to received as input:
#  - The path to the file to parse
#  - The information on the presence or not of column headers in the file
#  - The ordered list of column headers (if headers are not present or if user want to overwrite file header values)
#  - The name of the Class corresponding to Object to be created for each file line
#  - The ordered list of header names to use to call the Class constructor
#  - The ordered list of optional values to use to call the Class constructor
#  - The comment symbol
#  - The information indicating if the related db table must be cleaned before insertion or not
class FastaParser(object):

    ##
    # This method parses the provided TSV file, looking for headers if required and creating one Object of the requested Class per line
    #
    # @param file_path : string - The path to the TSV file to parse
    # @param regular_expression : string - The regular expression that must match the sequence definition line
    # @param class_name : string - name of the class to instantiate
    # @param parameter_name_list : list<string> - Ordered list of parameters used to instantiate the Class
    # @param parameter_value_list : list<string> - Ordered list of optional parameter values used to instantiate the Class
    # @param comment_symbol : string - symbol used to comment lines in the file to parse
    # @param clean_table : boolean  - The information indicating if the related db table must be cleaned before insertion or not
    #
    # @return None
    @staticmethod
    def parse_file( file_path, class_name, regular_expression, group_list, parameter_name_list, parameter_value_list, comment_symbol, clean_table = True):
        
        # Initialize the status of insertion
        status = Constants.STATUS_OK
        
        # Open the file to parse in 'read' mode
        input_file = FileUtils.open_text_r( file_path)                
        
        # Build the map of parameter name to regex group index
        # If a parameter name is not in the regex group, the index is set to -1
        # and a warning message is sent
        parameter_to_group_index_map = {}
        for parameter_name in parameter_name_list:
            try:
                index = group_list.index( parameter_name)
            except ValueError:
                index = -1
                Logger.get_instance().warning( "FastaParser.parse_file : The parameter '" + parameter_name + "' is not in the group list : " + str(group_list) + ".\n Check whether it is normal or not.")
                status = Constants.STATUS_WARNING
            parameter_to_group_index_map[ parameter_name] = index
        
        # Build the Factory to be used to create the objects
        data_factory = DataFactory( class_name)
        if clean_table:
            data_factory.clean_table()
        
        # compile the regular expression
        pattern = re.compile( regular_expression)
        
        # Parse the lines in the files and create one object by line matching the provided regex and
        # add it the sequence (lines below the matching line and before the next matching line or EOF)
        instance_list = []
        line = input_file.readline()
        current_matcher = None
        end_of_file = False
        while not end_of_file:
            # Read next line
            next_line = input_file.readline()
            if next_line == None or next_line == '':
                end_of_file = True
            # Ignore the empty and comment lines
            if (line != '' and not line.startswith( comment_symbol)) or end_of_file:
                # Check if the line match the sequence definition pattern
                matcher = pattern.match( line)
                if matcher != None or end_of_file:
                    # If a previous sequence definition already exists, create the corresponding object
                    if current_matcher != None or end_of_file:
                        # Build the ordered list of value required for Object creation
                        # If no parameter_value_list is provided, create a new 'None' one
                        # If a parameter_value_list is provided, override only the values
                        # corresponding to a parameter found in the file header
                        if( parameter_value_list == None):
                            parameter_value_list = [None]*len( parameter_name_list)
                        value_index = 0
                        for parameter_name in parameter_name_list:
                            parameter_index = parameter_to_group_index_map[ parameter_name]
                            if parameter_index >=0 :
                                parameter_value_list[ value_index] = current_matcher.group(group_list[ parameter_index])
                            value_index = value_index + 1
                        
                        #Call for the Object creation
                        new_instance = data_factory.create_object_from_fasta( parameter_value_list)
                        if new_instance != None:
                            instance_list.append( new_instance)
                    
                    # Initialize the new sequence definition and sequence
                    current_matcher = matcher
            
            line = next_line
            
        # Close the input file
        input_file.close()
        
        # Commit the SQLAlchemy session
        Logger.get_instance().info( "FastaParser.parse_file : Committing SQL session.")
        SQLManager.get_instance().commit();
        
        # Return the list of created instances
        return status
        
        
