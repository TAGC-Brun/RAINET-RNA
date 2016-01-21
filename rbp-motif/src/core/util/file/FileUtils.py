# -*- coding: utf-8 -*-

from core.util.log.Logger import Logger
from core.util.exception.FileFormatException import FileFormatException
from core.util.exception.RbpmotifException import RbpmotifException
from core.util import Constants

## This class contains static methods that are used to manipulate files like open
# and close, but some others like find extension.
# 
class FileUtils(object):


    #  open_text_r
    # --------
    #
    # Open a text file in reading mode.
    # Write a critical error in log file in case of IOError
    #
    # Argument(s):
    #    - path : the input file's path
    #
    # Return :
    #    - file_handle : is an object of type TODO
    @staticmethod
    def open_text_r( path):
        try:
            file_handle = open(path, 'r')
        except IOError as ioe:
            Logger.get_instance().critical\
            ("FileUtils.open_text_r : IOError:  Unable to open '" + path + "' : " + str(ioe))
            raise RbpmotifException( "Unable to open file '" + path + "' : " + str(ioe), ioe)

        return file_handle


    # open_text_w
    # -------
    #
    # Open a text file in writing mode
    # Write a critical error in log file in case of IOError
    #
    # Argument(s):
    #    - path : the input file's path
    #
    # Return :
    #    - file_handle : is an object of type TODO
    @staticmethod
    def open_text_w(path):
        try:
            file_handle = open(path, 'w')
        except IOError as detail:
            Logger.get_instance().critical\
            ("IOError:  Unable to open " + path + " : " + str(detail))
            exit()

        return file_handle
    
    
    
    # open_text_a
    # -------
    #
    # Open a text file in writing append mode
    # Write a critical error in log file in case of IOError
    #
    # Argument(s):
    #    - path : the input file's path
    #
    # Return :
    #    - file_handle : is an object of type TODO
    @staticmethod
    def open_text_a(path):
        try:
            file_handle = open(path, 'a')
        except IOError as detail:
            Logger.get_instance().critical\
            ("IOError:  Unable to open " + path + " : " + str(detail))
            exit()
        
        return file_handle
    
    
    

    # find_extension
    # -------------
    #
    # Find extension of the given file in path by searching in path's name
    # the extension : 'txt', 'ods', 'xls', 'xlsx'
    #
    # Return the extension in string uppercased format.
    #  
    @staticmethod
    def check_extension(path):
        # Split the path name on the '.'
        path_name_list = path.split('.')

        # Keep the extension ( last item on the list)
        extension = path_name_list[-1]

        # Transform all letters of extension in uppercase
        extension = extension.upper()

        if extension not in Constants.EXTENSION_LIST:
            raise FileFormatException\
            ("FileUTils.find_extension : FileFormatException : " + path + \
            " is not a file among with extension among " + Constants.EXTENSION_LIST +". Extension '" + extension + \
            "' is not recognized.")

        return extension
