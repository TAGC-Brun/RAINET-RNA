##
# Generic exception used to identify catched and raised exceptions and errors from Web code
#
class WebConnectionException( Exception):
    
    ##
    # Initialize the exception
    #
    # @param mess : string - The message associated to the exception
    # @param excep : Exception - The Exception associated to the exception (optional)
    def __init__(self, mess, excep = None):
        
        self.exception = excep
        self.message = mess

    ##
    # Display the exception message and stack trace
    def to_string(self):
        
        mess = self.__class__.__name__ + " :: " + self.message + " :\n"
        if self.exception != None:
            if self.exception.__class__.__name__ == "WebConnectionException":
                mess = mess + "  from : " + self.exception.to_string()
            else:
                mess = mess + str(self.exception)
                
        return mess
        