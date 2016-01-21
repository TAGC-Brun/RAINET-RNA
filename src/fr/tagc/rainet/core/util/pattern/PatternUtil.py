import re
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger



# #
# This class contains utilities about regular expression and pattern
class PatternUtil( object ):
    
    ##
    # Match the provided input string with the given compiled pattern
    # to return the group defined in the pattern.
    #
    # @param compiled_pattern : Pattern - A pattern compiled with the re.compile method. The regular expression
    #                            used to compile must contain a single group (without name)
    # @param input_string : string - The string to match with the Pattern
    #
    # @return the matched group if any, None if no match is found
    # @exception RainetException If more than one group is found in the matching
    @staticmethod
    def find_groups_in_string( compiled_pattern, input_string ):
        
        # Match the string to the pattern
        group_list = re.findall( compiled_pattern, input_string )
        
        # If there is at least a match, check for matched groups
        if group_list != None:
            return group_list

        return None
    
    ##
    # Retrieve the group in the query string provided a pattern with group definition
    #
    # @param pattern : string - The pattern to test. It must contain a group definition
    # @param query_string : string - The string to match on the pattern
    #
    # @return string - The group string if found or if not, either None if strict mode is true
    #                   or the query_string if strict mode is false
    #
    @staticmethod
    def find_single_group_in_string( pattern, query_string, strict = True):
        
        #Compile the pattern
        compiled_pattern = re.compile( pattern )
        # look for the group in the query_string matching the pattern
        group_list = PatternUtil.find_groups_in_string( compiled_pattern, query_string )
        # If there is group(s)
        if group_list != None:
            # If there is one group, return it
            if len( group_list) == 1:
                return group_list[0]
            # If there is more than group, warn user. 
            # If strict mode is true, return None, if it is false, return the query_string
            elif len( group_list) > 1 :
                Logger.get_instance().warning( "PatternUtil.find_single_group_in_string : string contains several groups : pattern = " + pattern + " and groups = " + str( group_list))
                if strict:
                    return None
                else:
                    return group_list
            # If there is no group, warn user.
            # If strict mode is true, return None, if it is false, return the query_string
            else:
                Logger.get_instance().warning( "PatternUtil.find_single_group_in_string : string contains no group : pattern = " + pattern + " and query string = " + query_string)
                if strict:
                    return None
                else:
                    return query_string
        
        # If there is no group, warn user.
        # If strict mode is true, return None, if it is false, return the query_string
        else:
            Logger.get_instance().warning( "PatternUtil.find_single_group_in_string : string contains no group : pattern = " + pattern + " and query string = " + query_string)
            if strict:
                return None
            else:
                return query_string
    
    ##
    # Match the strings of the provided input list with the given compiled pattern
    # to return the groups defined in the pattern. 
    #
    # @param compiled_pattern : Pattern - A pattern compiled with the re.compile method. The regular expression
    #                            used to compile must contain a single group (without name)
    # @param input_list : list<string> - The list of strings to be matched with the Pattern
    #
    # @return the list of groups found. May be empty.
    # @exception RainetException If no group is found and an exception was raised during the search
    @staticmethod
    def find_groups_in_list( compiled_pattern, input_list):
        
        # Parse the list of string to match
        exception_raised = None
        group_list = []
        for input_string in input_list:
            # Search for the group in the current string
            try:
                Logger.get_instance().debug( "find_group_in_list = " + input_string)
                groups = PatternUtil.find_groups_in_string( compiled_pattern, input_string)
                if groups != None and len( groups) > 0:
                    group_list.extend( groups)
            except RainetException as e:
                exception_raised = e
                    
        # If no group was found and an exception was raised during the process
        # raise it 
        if len(group_list) == 0 and exception_raised != None:
            raise exception_raised
        
        # If a group was found, return it
        return group_list
            