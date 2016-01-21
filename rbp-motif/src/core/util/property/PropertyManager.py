from ConfigParser import ConfigParser

from core.util.log.Logger import Logger
from core.util.exception.RbpmotifException import RbpmotifException







# # This class is a singleton aiming to manage the properties retrieved from properties files
# like insertion.ini file that contains list of files used for data MakeDatasetRbp.
#
class PropertyManager( object ) :

    __instance = None
    
    # #
    #
    def __init__( self ):
        
        self.propertiesDict = {}
        
    
    ##
    # Read the properties of a file and store them into a dictionary, removing the section information
    #
    # @param file_path : string - The path to the properties file
    def read_properties(self, file_path ):
        
        Logger.get_instance().info("Reading properties from file : " + file_path)
        config_parser = ConfigParser()
        
        config_parser.read( file_path)

        for section in config_parser.sections():
            options = config_parser.options(section)
            for option in options:
                
                try:
                    option = option.lower()
                    self.propertiesDict[option] = config_parser.get(section, option)
                    if self.propertiesDict[option] == None:
                        raise RbpmotifException( "PropertyManager.readProperties : Wrong property definition for key = " + option + " in section " + section + " of file " + file_path)
                except:
                    raise RbpmotifException( "PropertyManager.readProperties : Abnormal property definition for key = " + option + " in section " + section + " of file " + file_path)

    # #
    # Returns the value of the property with the provided property_name.
    # if the property is not available, return None or an Exception
    #
    # @param property_name : string - The name of the property to get
    # @param not_none : boolean - (optional) If True, a None value can be returned when
    #                     property is not available. If False, an exception is raised in such case
    # @return The value of the property if present.
    # @raise RainetException : When the property is not available and not_none is set to True 
    def get_property(self, property_name, not_none = False):
        
        property_name_lower = property_name.lower()
        if property_name_lower in self.propertiesDict.keys():
            return self.propertiesDict[property_name_lower]
        else:
            if not_none:
                print str(self.propertiesDict)
                raise RbpmotifException( "PropertyManager.get_property : No property with name : " + property_name_lower)
            else:
                return None

    # #
    # Returns the singleton instance
    #
    # @return the singleton instance
    @staticmethod
    def get_instance():

        if PropertyManager.__instance == None:
            PropertyManager.__instance = PropertyManager()
        return PropertyManager.__instance