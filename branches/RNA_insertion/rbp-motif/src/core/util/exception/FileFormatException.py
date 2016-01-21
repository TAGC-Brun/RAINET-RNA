# -*- coding: utf-8 -*-




##
# Exception raised in case of file format given in a method is not
# appropriated.
#
class FileFormatException():


    #  Constructor of FileFormatException
    # -----------------------------------
    #
    # Arguments :
    #    - message : a message describing in which context the exception happened.
    def __init__ (self, message):
        self.message = message

    # print_error
    # ----------
    #
    #  Return a string explaining the exception.
    # The string is the message given when the exception is raised.
    def print_error(self):
        string = self.message
        return string
