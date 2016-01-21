
import re

from sqlalchemy import or_

from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data.NetworkModule import NetworkModule
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.util import Constants
from fr.tagc.rainet.core.util.pattern.PatternUtil import PatternUtil
from fr.tagc.rainet.core.data.NetworkModuleAnnotation import NetworkModuleAnnotation
from fr.tagc.rainet.core.data.ProteinNetworkModule import ProteinNetworkModule
from fr.tagc.rainet.core.util.factory.DataFactory import DataFactory
from fr.tagc.rainet.core.data import DataConstants

## This class is a parser of a graph module (class) annotation file produced by the Brun team tools
# Each retrieved module get its annotations added
# The parser have to received as input:
# -file_path : string - The fpath to the file to parse
# -class_tag : string - The tag indicating the beginning of a class (module) definition
# -class_regex : string - The regex that provides the location of the class (module) ID as a group
# -protein_tag : string - The tag indicating the list of class (module) proteins
# -annotation_tag : string - The tag indicating the list of class (module) annotations
# -comment_char : string - The character used to declare comment lines
# -clean_table : boolean  - The information indicating if the related db table must be cleaned before insertion or not
class NetworkModuleAnnotationParser ( object ):

    # #
    # Parse the network module file and build the NetworkModule objects. Clean the NetworkModuleAnnotation table if required
    #
    # @param file_path : string - The fpath to the file to parse
    # @param class_tag : string - The tag indicating the beginning of a class (module) definition
    # @param class_regex : string - The regex that provides the location of the class (module) ID as a group
    # @param protein_tag : string - The tag indicating the list of class (module) proteins
    # @param annotation_tag : string - The tag indicating the list of class (module) annotations
    # @param comment_char : string - The character used to declare comment lines
    # @param clean_table : boolean  - The information indicating if the related db table must be cleaned before insertion or not
    #
    # @return Constants.STATUS_OK if everything went fine
    # @raise RainetException if an error occured during parsing
    @staticmethod
    def parse_file( file_path, class_tag, class_regex, protein_tag, annotation_tag, comment_char, clean_table = True ):
        
        # Build the Factory to Clean the table
        data_factory = DataFactory( DataConstants.INTERACTOME_NETWORK_PARTITION_ANNOTATION_CLASS)
        if clean_table:
            data_factory.clean_table()
        
        # Parse the file to get the modules
        NetworkModuleAnnotationParser.process_modules( file_path, class_tag, class_regex, protein_tag, annotation_tag, comment_char )
       
        # Commit the SQLAlchemy session
        Logger.get_instance().info( "NetworkModuleAnnotationParser.parse_file : Committing SQL session." )
        SQLManager.get_instance().commit();
        
        return Constants.STATUS_OK
        

    # #
    # Parses a Network module annotation file
    #
    # @param input_file : string - The fpath to the file to parse
    # @param class_tag : string - The tag indicating the beginning of a class (module) definition
    # @param class_regex : string - The regex that provides the location of the class (module) ID as a group
    # @param protein_tag : string - The tag indicating the list of class (module) proteins
    # @param annotation_tag : string - The tag indicating the list of class (module) annotations
    # @param comment_char : string - The character used to declare comment lines
    #
    # @return None
    # @raise RainetException if an error occurred during parsing
    @staticmethod
    def process_modules( input_file, class_tag, class_regex, protein_tag, annotation_tag, comment_char ):

        sql_session = SQLManager.get_instance().get_session()
        
        class_pattern = re.compile( class_regex )
        
        with open( input_file, "r" ) as infile:
            current_module = None
            for line in infile:
                # TODO remove print
                line = line.strip()
                # Skip empty line
                if not line:
                    continue
                # Skip comment line
                if line[0] == comment_char:
                    continue
                # Module (class) line
                elif line.startswith( class_tag):
                    # If there is a previous defined NetworkModule, insert it in DB
                    if current_module != None:
                        sql_session.add( current_module )
                    # Start the new NetworkModule
                    module_id = PatternUtil.find_single_group_in_string( class_pattern, line )
                    # -- Find the corresponding module
                    current_module_list = sql_session.query( NetworkModule).filter( NetworkModule.moduleID == module_id).all()
                    if current_module_list == None or len( current_module_list) == 0:
                        raise RainetException( "NetworkModuleAnnotationParser.process_modules : a module ID was not found in DB : " + line)
                    if len( current_module_list) > 1:
                        raise RainetException( "NetworkModuleAnnotationParser.process_modules : several NetworkModuel found in DB for module ID " + module_id + " (file line =" + line + ")")
                    current_module = current_module_list[0]
                # Protein list associated to module
                elif line.startswith( protein_tag): 
                    if current_module != None:
                        NetworkModuleAnnotationParser.check_proteins( current_module, line, sql_session )
                    else:
                        raise RainetException( "NetworkModuleAnnotationParser.process_modules : Trying to associate protein while Network module is None : " + line )
                # Annotation list associated to module
                elif line.startswith( annotation_tag): 
                    if current_module != None:
                        NetworkModuleAnnotationParser.associate_annotations( current_module, line, sql_session )
                    else:
                        raise RainetException( "NetworkModuleAnnotationParser.process_modules : Trying to associate protein while Network module is None : " + line )                    
            # Add last Module
            if current_module != None:
                sql_session.add( current_module )


    ##
    # Check if the list of proteins provided in the list correspond to the network module proteins listed in the DB
    #
    # @param network_module : NetworkModule - The network module the proteins should be part of
    # @param line : string - the line of the file containing the list of proteins
    # @param sql_session: the current SQLAlchemy session
    #
    # @return None
    # @raise RainetException if an error occurred during parsing
    @staticmethod
    def check_proteins( network_module, line, sql_session ):
        
        # Get all the proteins AC associated to the networkModule
        # to compare them with the proteins in the annotated module
        module_proteins = network_module.associatedProteins
        module_protein_id_list = [] 
        for protein in module_proteins:
            module_protein_id_list.append( protein.uniprotAC)
        
        # Split the line with tab to separate tag and protein list 
        token_list = line.split( "\t")
        if token_list == None or len( token_list) != 2:
            raise RainetException( "NetworkModuleAnnotationParser.check_proteins : Abnormal annotation line : " + line + " in module " + network_module.moduleID)
        
        # Get the protein list by splitting with spaces
        protein_list = token_list[1].split(", ")
        
        # Parse the proteins to check if they are associated to the given network module    
        for protein_name in protein_list:
            # Get the protein with the protein_name as ID (uniprotAD or cross reference)
            # Search first in Protein table both with Uniprot AC and Uniprot ID
            # then with CrossReference
            protein_list = sql_session.query( Protein ).filter( or_( Protein.uniprotID == protein_name, Protein.uniprotAC == protein_name ) ).all()
            if protein_list == None or len( protein_list ) == 0:
                protein_list = sql_session.query( Protein ).filter( Protein.uniprotAC == ProteinCrossReference.protein_id, ProteinCrossReference.crossReferenceID == protein_name ).all()
            
            # Check if a single Protein is found. If not, trace issue and continue.
            # If so remove the protein from the list of proteins found in the module
            # to account for consistency
            protein = None
            if protein_list != None and len( protein_list ) > 0:
                if len( protein_list ) == 1:
                    protein = protein_list[0]
                    if protein != None:
                        if protein.uniprotAC in module_protein_id_list:
                            module_protein_id_list.remove( protein.uniprotAC)
                        else:
                            raise RainetException( "NetworkModuleAnnotationParser.check_proteins : Consistency error : a protein of the annotated module '" + network_module.moduleID + "' was not found in the Network module : " + protein.uniprotAC + "/" + protein_name)
                    
                else:
                    Logger.get_instance().debug( "NetworkModuleAnnotationParser.check_proteins : Abnormal number of Protein found for ID = " + protein_name + " : " + str( len( protein_list ) )  + " in network module ID = " + network_module.moduleID )
                    continue
                
            # If no protein was found for the corresponding name, continue
            if protein == None:
                Logger.get_instance().debug( "NetworkModuleAnnotationParser.check_proteins : No Protein found for ID = " + protein_name + " in network module ID = " + network_module.moduleID )
                continue
            
            # Search for the relation from the protein to the network module
            #protein_network_module_list = sql_session.query( ProteinNetworkModule).filter( ProteinNetworkModule.protein_id == protein.uniprotAC, ProteinNetworkModule.networkModule_id == network_module.moduleID).all()
            #if protein_network_module_list == None or len( protein_network_module_list) == 0:
            #    raise RainetException( "NetworkModuleAnnotationParser.check_proteins : No relationship found in DB between protein " + protein.uniprotAC + " and network module ID " + network_module.moduleID) 
         
        # If the proteins of the network module were not all found in the annotated module, raise an issue
        if len( module_protein_id_list) > 0:
            raise RainetException( "NetworkModuleAnnotationParser.check_proteins : Consistency error : some proteins of the Network module '" + network_module.moduleID + "' were not found in the annotated module : " + str( module_protein_id_list))
            
            
    ##
    # Parse the provided line to find annotations. Create corresponding NetworkModuleAnnotation objects
    # and attach them to the provided network_module
    #
    # @param network_module : NetworkModule - The annotated network module
    # @param line : string - the line of the file containing the list of annotations
    # @param sql_session: the current SQLAlchemy session
    #
    # @return None
    # @raise RainetException if an error occurred during parsing
    @staticmethod
    def associate_annotations( network_module, line, sql_session):
        
        # Split the line with tab to separate tag and annotation list 
        token_list = line.split( "\t")
        if token_list == None or len( token_list) != 2:
            raise RainetException( "NetworkModuleParser.associate_annotations : Abnormal annotation line : " + line + " for network module ID = " + network_module.moduleID)
        
        # Get the annotation list by splitting with spaces
        annotation_list = token_list[1].split(" ")
        
        # Parse the annotations and create the new NetworkModuleAnnotation object
        # adding it to the annotated network module    
        for annotation in annotation_list:
            current_annotation = NetworkModuleAnnotation( annotation)
            network_module.add_associated_annotation( current_annotation)
            sql_session.add( current_annotation)
            sql_session.add( network_module)
            
            
            
            
