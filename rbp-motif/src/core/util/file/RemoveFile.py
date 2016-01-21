# This class removes a file through os python module 


import os 
from core.util.log.Logger import Logger

class RemoveFile():
    
    # This method Removes a file
    # Argument:
    #    - file path
    @staticmethod
    def delete_file(namefile):
        
        try:
            os.remove(namefile)
            Logger.get_instance().info(" This file has been removed: " + namefile)
        except OSError:
            Logger.get_instance().info(" Cannot remove : " + namefile)
        
        


