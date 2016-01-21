
import os
from sqlalchemy import or_

from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data.NetworkModule import NetworkModule
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.PPINetwork import PPINetwork
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.PartitionAnalysis import PartitionAnalysis
from fr.tagc.rainet.core.util import Constants
from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.factory.DataFactory import DataFactory
from fr.tagc.rainet.core.util.file.FileUtils import FileUtils

## This class is a parser of a graph module (class) file produced by the Brun team tools
# Each retrieved module produces a NetworkModule and geta ssociated to the Protein that are part of it
# The parser have to received as input:
# - file_path : string - The fpath to the file to parse
# - class_tag : string - The tag indicating the beginning of a class (module) definition
# - comment_char : string - The character used to declare comment lines
# - clean_table : boolean  - The information indicating if the related db table must be cleaned before insertion or not
class NetworkModuleParser ( object ):

    # #
    # Parse the network module file and build the NetworkModule objects. Clean the NetworkModule table if required.
    #
    # @param file_path : string - The fpath to the file to parse
    # @param class_tag : string - The tag indicating the beginning of a class (module) definition
    # @param comment_char : string - The character used to declare comment lines
    # @param clean_table : boolean  - The information indicating if the related db table must be cleaned before insertion or not
    @staticmethod
    def parse_file( file_path, class_tag, comment_char, clean_table = True ):
        
        # Build the Factory to Clean the table
        data_factory = DataFactory( DataConstants.INTERACTOME_NETWORK_PARTITION_CLASS)
        if clean_table:
            data_factory.clean_table()
        
        # Parse the file to get the modules
        NetworkModuleParser.process_modules( file_path, class_tag, comment_char )
       
        # Commit the SQLAlchemy session
        Logger.get_instance().info( "NetworkModuleParser.parse_file : Committing SQL session." )
        SQLManager.get_instance().commit();
        
        return Constants.STATUS_OK
        

    # #
    # Parses a Network module file (.clas file)
    #
    # @param file_path : string - The fpath to the file to parse
    # @param class_tag : string - The tag indicating the beginning of a class (module) definition
    # @param comment_char : string - The character used to declare comment lines
    @staticmethod
    def process_modules( input_file, class_tag, comment_char ):

        sql_session = SQLManager.get_instance().get_session()
        
        # Create a new encapsulating OCGParittionAnalysis
        partition_analysis = PartitionAnalysis( os.path.basename( input_file ), "" )
        
        # Look for a PPINetwork that seems to correspond to the class file
        class_file_name = os.path.basename( input_file)
        no_extension_name = FileUtils.remove_extension( class_file_name )
        gr_extension_name = no_extension_name + ".gr"
        query = sql_session.query( PPINetwork).filter( or_( PPINetwork.name == class_file_name,
                                                            PPINetwork.name == no_extension_name,
                                                            PPINetwork.name == gr_extension_name)).all()
        
        if( query == None or len( query) == 0):
            raise RainetException( "NetworkModuleParser:process_modules : Unable to find a PPINetwork using name similar to " + class_file_name + " using " + no_extension_name + " or " + gr_extension_name)
        
        if( len( query) > 1):
            raise RainetException( "NetworkModuleParser:process_modules : Several (" + str( len( query)) +") PPINetwork found using name similar to " + class_file_name)
        
        ppi_network = query[0]
                
        # Add the PartitionAnalysis to PPINetwork
        ppi_network.add_partition_analysis( partition_analysis)
            
        # Parse the file
        comment = ""
        with open( input_file, "r" ) as infile:
            current_module = None
            for line in infile:
                line = line.strip()
                # Skip empty lines
                if not line: 
                    continue
                # Comment line
                if line[0] == comment_char:
                    comment = comment + line
                # Module (class) line
                elif line[0:len( class_tag)] == class_tag:
                    # If there is a previous defined NetworkModule, insert it in DB
                    if current_module != None:
                        sql_session.add( current_module )
                        Logger.get_instance().info( "Inserting new NetworkModule : " + current_module.moduleID)
                    # Start the new NetworkModule
                    # -- Locate the parenthesis after the class ID (if any)
                    parenthesis_index = line.index( "(", 6 )
                    # -- If there is a parenthesis, the ID is before
                    if parenthesis_index > 0:
                        module_id = line[ 7 : ( parenthesis_index - 1 )].strip()
                    # -- If there is no parenthesis, the ID is the final part of the line
                    else:
                        module_id = line[7:].strip()
                    # -- Instantiate the new module and add it to the partition analysis
                    current_module = NetworkModule( module_id )
                    partition_analysis.add_network_module( current_module)
                # Line not empty, not comment, not class : it is a module protein list or a line to ignore
                else: 
                    if current_module != None:
                        NetworkModuleParser.associate_proteins( current_module, line )
                    else:
                        Logger.get_instance().debug( "NetworkModuleParser.process_modules : Ignoring line = " + line )
                    
            # Add last Module to the SQL Session
            if current_module != None:
                sql_session.add( current_module )
        
        # Add the class file comment as analysis details
        partition_analysis.analysisDetails = comment
        
        # Insert PartitionAnalysis and PPINetwork to SQL session 
        sql_session.add( partition_analysis)
        sql_session.add( ppi_network)

    ##
    # Parse the provided line to find protein names. Create corresponding association between the proteins and the module
    #
    # @param module : NetworkModule - The network module the proteins must be associated to
    # @param line : string - the line of the file containing the list of protein names
    #
    # @return None
    # @raise RainetException if an error occurred during parsing
    @staticmethod
    def associate_proteins( module, line ):
        
        token_list = line.split( " " )
        
        sql_session = SQLManager.get_instance().get_session()
        
        if token_list != None and len( token_list ) > 0:
            for token in token_list:
                # Get the protein with the token as ID (uniprotAD or cross reference)
                # Search first in Protein table both with Uniprot AC and Uniprot ID
                # then with CrossReference
                protein_list = sql_session.query( Protein ).filter( or_( Protein.uniprotID == token, Protein.uniprotAC == token ) ).all()
                if protein_list == None or len( protein_list ) == 0:
                    protein_list = sql_session.query( Protein ).filter( Protein.uniprotAC == ProteinCrossReference.protein_id, ProteinCrossReference.crossReferenceID == token ).all()
                # Check if a single Protein is found. If not, raise an issue or report problem
                protein = None
                if protein_list != None and len( protein_list ) > 0:
                    if len( protein_list ) == 1:
                        protein = protein_list[0]
                        sql_session.add( protein )
                    else:
                        raise RainetException( "NetworkModuleParser.associate_proteins : Abnormal number of Protein found for ID = " + token + " : " + str( len( protein_list ) ) ) 
                else:
                    Logger.get_instance().debug( "NetworkModuleParser.associate_proteins : No Protein found for ID = " + token )
                    continue
                # if a protein is found, associate it to the module
                if protein != None:
                    module.add_associated_protein( protein )
                else:
                    raise RainetException( "NetworkModuleParser.associate_proteins : Protein found is none for ID = " + token )
        else:
            raise RainetException( "NetworkModuleParser.associate_proteins : list of proteins is empty or None : " + line )
        
