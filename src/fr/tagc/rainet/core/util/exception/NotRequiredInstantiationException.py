
import traceback

from fr.tagc.rainet.core.util.exception.RainetException import RainetException

##
# This Exception is used when automatic instantiation is used through Factory and
# the class instantiation is not required since a similar object already exist (in DB for instance)
#
class NotRequiredInstantiationException( RainetException):
    
    ##
    # Display the exception message and stack trace
    def to_string(self):
        
        mess = "NotRequiredInstantiationException : " + self.message + " : " + str( self.exception) + "\n" + traceback.format_exc()
        return mess