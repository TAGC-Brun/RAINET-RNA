#  ***** Class access web ******




from core.util.log.Logger import Logger
from urllib2 import urlopen
from core.util.exception.WebConnectionException import WebConnectionException
from core.util import Constants
import time

class Web():    
    
    # openurl
    # -------------------------
    #     
    # This function access to specific site in a monitored way
    # Arguments:
    #    - url
    #    - text = True --> returns the url page in string format
    #    - text = False --> returns the url page in page format
    #
    @staticmethod
    def openurl(url, text=True):
        connected = False
        count = 0
        while connected == False and count <= Constants.WEB_CONNECTION_CONSTANT:
            try:
                count += 1
                time.sleep( 1)
                print "Trying to connect : count =" + str( count)
                page = urlopen(url, timeout=3)
                connected = True
                print 'good connection'
                if text==True:
                    print "Get text"
                    string = ''
                    for line in page:
                        string+=line
                    print "Return line"
                    return string
                else:
                    print "Return page"
                    return page
                time.sleep( count)
            except IOError: #as ioe: # Aggiungere altre exception
                print 'BAD connection'
                Logger.get_instance().warning(" Web.openurl : Problem to establish a connection to " + url) 
                connected = False
                #raise WebConnectionException( "Web.openurl : Problem to establish a connection to " + url, ioe)
            except Exception as e:
                print 'BAD connection'
                Logger.get_instance().warning(" Web.openurl : Unexpected Error : " + e)
                connected = False
        
        print "End Try"
        
        if connected == False:
            Logger.get_instance().error(" Web.openurl : unable to connect to url" + url + " after "+ str(Constants.WEB_CONNECTION_CONSTANT) + " retries" )
            raise WebConnectionException( "Web.openurl : unable to establish a connection to " + url)
            
        