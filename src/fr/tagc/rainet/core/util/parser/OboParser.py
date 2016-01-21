
from __future__ import with_statement
from collections import defaultdict

from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data.GeneOntology import GeneOntology
from fr.tagc.rainet.core.util import Constants
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.factory.DataFactory import DataFactory
from fr.tagc.rainet.core.data import DataConstants

__author__ = "Uli Koehler"
__copyright__ = "Copyright 2013 Uli Koehler"
__license__ = "Apache v2.0"

## This class is a parser of the Gene Ontology definition file.
# Each retrieved GO Term produces a GeneOntology object
# The parser have to received as input:
# - file_path : string - The fpath to the file to parse
# - go_id_tag : string - The tag declaring the GO term ID
# - co_name_tag : string - The tag declaring the GO term name
# - go_namespace_tag : string - The tag declaring the GO term namespace
# - clean_table : boolean  - The information indicating if the related db table must be cleaned before insertion or not
class OboParser ( object ):


    ##
    # Parse the obo file and build the GeneOntology objects
    #
    # @param file_path : string - The path of the file to parse
    # @param go_id_tag : string - The tag declaring the GO term ID
    # @param co_name_tag : string - The tag declaring the GO term name
    # @param go_namespace_tag : string - The tag declaring the GO term namespace
    # @param clean_table : boolean - The information indicating if the related db table must be cleaned before insertion or not
    #
    @staticmethod
    def parse_file( file_path, go_id_tag, go_name_tag, go_namespace_tag, clean_table = True):
        
        # Parse the file
        go_term_list = OboParser.processGOOBO( file_path)

        # Build the Factory to clean the table
        data_factory = DataFactory( DataConstants.GENE_ONTOLOGY_CLASS)
        if clean_table:
            data_factory.clean_table()
        
        # Get the SQL session
        sql_session = SQLManager.get_instance().get_session()
        
        # Build the GeneOntology objects from the parsing result
        for go_term in go_term_list:
            go_id = go_term[ go_id_tag]
            go_name = go_term[ go_name_tag]
            go_namespace = go_term[ go_namespace_tag]
            if go_id != None and go_id != '' and go_name != None and go_name != '' and go_namespace != None and go_namespace != '':
                Logger.get_instance().debug( "OboParser.parse_file : inserting GO " + go_id + " " + go_name + " " + go_namespace)
                new_go = GeneOntology( go_id, go_name, go_namespace)
                sql_session.add( new_go)
                
       
        # Commit the SQLAlchemy session
        Logger.get_instance().info( "OboParser.parse_file : Committing SQL session.")
        SQLManager.get_instance().commit();
        
        return Constants.STATUS_OK
        

    ##
    #In an object representing a GO term, replace single-element lists with
    #their only member.
    # 
    # @param goTerm : string - The GO Term
    # @return The modified object as a dictionary
    @staticmethod
    def processGOTerm( goTerm ):
        ret = dict( goTerm )  # Input is a defaultdict, might express unexpected behaviour
        for key, value in ret.iteritems():
            if len( value ) == 1:
                ret[key] = value[0]
        return ret

    ##
    # Parses a Gene Ontology dump in OBO v1.2 format.
    # Yields each 
    #
    # @param input_file: string - The filename to read
    @staticmethod
    def processGOOBO( input_file ):

        with open( input_file, "r" ) as infile:
            currentGOTerm = None
            for line in infile:
                line = line.strip()
                if not line: continue  # Skip empty
                if line == "[Term]":
                    if currentGOTerm: yield OboParser.processGOTerm( currentGOTerm )
                    currentGOTerm = defaultdict( list )
                elif line == "[Typedef]":
                    # Skip [Typedef sections]
                    currentGOTerm = None
                else:  # Not [Term]
                    # Only process if we're inside a [Term] environment
                    if currentGOTerm is None: continue
                    key, sep, val = line.partition( ":" )
                    currentGOTerm[key].append( val.strip() )
            # Add last term
            if currentGOTerm is not None:
                yield OboParser.processGOTerm( currentGOTerm )
  
