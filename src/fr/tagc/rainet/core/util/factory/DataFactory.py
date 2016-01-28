
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException

from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinDomain import ProteinDomain
from fr.tagc.rainet.core.data.ProteinInteraction import ProteinInteraction
from fr.tagc.rainet.core.data.ProteinIsoform import ProteinIsoform
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.KEGGPathway import KEGGPathway
from fr.tagc.rainet.core.data.ProteinKEGGAnnotation import ProteinKEGGAnnotation
from fr.tagc.rainet.core.data.ReactomePathway import ReactomePathway
from fr.tagc.rainet.core.data.ProteinReactomeAnnotation import ProteinReactomeAnnotation
from fr.tagc.rainet.core.data.GeneOntology import GeneOntology
from fr.tagc.rainet.core.data.ProteinGOAnnotation import ProteinGOAnnotation
from fr.tagc.rainet.core.data.NetworkModule import NetworkModule
from fr.tagc.rainet.core.data.NetworkModuleAnnotation import NetworkModuleAnnotation
from fr.tagc.rainet.core.data.PPINetwork import PPINetwork
from fr.tagc.rainet.core.data.PPINetworkInteraction import PPINetworkInteraction
from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.RNACrossReference import RNACrossReference


## This class is a Factory that aims to create objects of various Class from class name
# and required parameters. The class manage the insertion to database of the newly created
# object and ensure its correct connection in the database schema.
class DataFactory:

    # #
    # Constructor of the Factory
    #
    # @param class_name : string - The name of the Class to instantiate 
    #
    def __init__( self, class_name ):

        self.className = class_name

    ##
    # Remove the content of the given table
    #
    def clean_table(self):
        
        Logger.get_instance().info( "|--Cleaning table : " + self.className)
        sql_session = SQLManager.get_instance().get_session()
        line_list = eval( "sql_session.query(" + self.className + ").all()" )
        if line_list != None:
            for line in line_list:
                sql_session.delete( line)
            SQLManager.get_instance().commit()
            Logger.get_instance().info( "|--Table cleaned : " + self.className)
        else:
            Logger.get_instance().warning( "|--Table not found : " + self.className)                
        
    ##
    # Create an instance of the given class with the given parameters using the correct constructor
    # Insert the instance to database using the SQL session provided at Factory creation
    #
    # @param class_name : string - the name of the Class to instantiate
    # @param parameter_value_list : list - the list of parameters values to be used to create the object
    # 
    # @raise RainetException if an error occurred during object creation.
    #
    # @return The created instance
    #
    # @raise RainetException if an error occurred while building or inserting the new instance to the database
    def create_object_from_tsv( self, parameter_value_list ):
        # Test if the parameter list contains elements
        if parameter_value_list == None or len( parameter_value_list ) == 0:
            raise RainetException( "DataFactory.create_object_from_tsv: No value provided for the object creation of class " + self.className )
        
        #=======================================================================
        # Build the code (introspection) for constructor call
        #=======================================================================
        # -- start with the variable that will receive the new instance
#         constructor_command = "new_instance = "
        constructor_command = ""
        # -- add the class name with one parenthesis
        constructor_command += self.className + "("
        # -- add the list of parameters separated by a comma
        for parameter_value in parameter_value_list:
            constructor_command += "\"" + str( parameter_value.replace("\"", "\\\"")) + "\","
        # -- remove the last comma
        constructor_command = constructor_command[:-1]
        # -- add the final closing parenthesis
        constructor_command += ")"
        
        #=======================================================================
        # Evaluate the constructor command and add the object to SQLAlchemy session if required
        #=======================================================================
        new_instance = None
        try:
            Logger.get_instance().debug( "command = " + constructor_command)
            new_instance = eval( constructor_command )
            new_instance.add_to_session()
        # Possibly the object does not have to be created because it is already in database
        # In that case, a NotRequiredInstantiationException is raised by the class constructor
        except NotRequiredInstantiationException:
            pass
        # If any other Exception occurred, something went wrong
        except Exception as excep:
            raise RainetException( "DataFactory.create_object_from_tsv : An issue occurred in the object creation.", excep )
        
        return new_instance

    # #
    # Create an instance of the given class with the given parameters using the correct constructor
    # Insert the instance to database using the SQL session provided at Factory creation
    #
    # @param class_name : string - the name of the Class to instantiate
    # @param parameter_value_list : list - the list of parameters values to be used to create the object
    # 
    # @raise RainetException if an error occurred during object creation.
    #
    # @return The created instance
    #
    # @raise RainetException if an error occurred while building or inserting the new instance to the database
    def create_object_from_fasta( self, parameter_value_list):
        # Test if the parameter list contains elements
        if parameter_value_list == None or len( parameter_value_list ) == 0:
            raise RainetException( "DataFactory.create_object_from_fasta: No value provided for the object creation of class " + self.className )
        
        #=======================================================================
        # Build the code (introspection) for constructor call
        #=======================================================================
        # -- start with the variable that will receive the new instance
        # constructor_command = "new_instance = "
        constructor_command = ""
        # -- add the class name with one parenthesis
        constructor_command += self.className + "("
        # -- add the list of parameters separated by a comma
        for parameter_value in parameter_value_list:
            constructor_command += "\"" + str(parameter_value) + "\","
        # -- remove the last comma
        constructor_command = constructor_command[:-1]
        # -- add the final closing parenthesis
        constructor_command += ")"
        
        #=======================================================================
        # Evaluate the constructor command and add the object to SQLAlchemy session if required
        #=======================================================================
        new_instance = None
        try:
            Logger.get_instance().debug( "command = " + constructor_command)
            new_instance = eval( constructor_command )
            new_instance.add_to_session()
        # Possibly the object does not have to be created because it is already in database
        # In that case, a NotRequiredInstantiationException is raised by the class constructor
        except NotRequiredInstantiationException:
            pass
        # If any other Exception occurred, something went wrong
        except Exception as excep:
            raise RainetException( "DataFactory.create_object_from_fasta : An issue occurred in the object creation.", excep )
        
        return new_instance