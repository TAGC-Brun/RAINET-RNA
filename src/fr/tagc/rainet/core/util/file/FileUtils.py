# -*- coding: utf-8 -*-

import os

from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.exception.FileFormatException import FileFormatException
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util import Constants
from string import rfind

## This class contains static methods that are used to manipulate files like open
# and close, but some others like find extension.
# 
class FileUtils(object):


    ## initialise_output_folder
    #  --------------
    #
    # Initialise all Analysis output folders
    #
    # Return nothing
    @staticmethod
    def initialise_output_folders( base_folder):

        # If base folder does not exist
        if not os.path.exists( base_folder):
            os.mkdir( base_folder)

        # Try to create report folder
        if not os.path.exists( base_folder+"/"+Constants.REPORT_FOLDER):
            os.mkdir( base_folder+"/"+Constants.REPORT_FOLDER)

        return


    #  open_text_r
    # --------
    #
    # Open a text file in reading mode.
    # Write a critical error in log file in case of IOError
    #
    # @param path : the input file's path
    #
    # @return an object of type file handler
    @staticmethod
    def open_text_r(path):
        try:
            file_handle = open(path, 'r')
        except IOError as ioe:
            Logger.get_instance().critical\
            ("FileUtils.open_text_r : IOError:  Unable to open '" + path + "' : " + str(ioe))
            raise RainetException( "Unable to open file '" + path + "' : " + str(ioe), ioe)

        return file_handle


    # open_text_w
    # -------
    #
    # Open a text file in writing mode
    # Write a critical error in log file in case of IOError
    #
    # @param path : string - the input file's path
    #
    # @return :
    #    - file_handle : is an object of type TODO
    @staticmethod
    def open_text_w(path):
        try:
            file_handle = open(path, 'w')
        except IOError as detail:
            Logger.get_instance().critical\
            ("IOError:  Unable to open " + path + " : " + str(detail), detail)
            exit()

        return file_handle


    # find_extension
    # -------------
    #
    # Find extension of the given file in path by searching in path's name
    # the extension : 'txt', 'ods', 'xls', 'xlsx'
    #
    # @param path : string - the file path
    #
    # @return the extension in string uppercase format.
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
    
    ##
    # Remove the file extension from the file name if it exists
    #
    # @param file_name : string - the file name or file path
    #
    # @return the file_name without the file extension
    @staticmethod
    def remove_extension( file_name):
        
        dot_index = rfind( file_name, ".")
        
        if( dot_index > 0):
            return file_name[0:dot_index]
        else:
            return file_name
