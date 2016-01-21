
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.TableStatus import TableStatus
from fr.tagc.rainet.core.util import Constants
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from sqlalchemy.sql.schema import MetaData
from sqlalchemy import or_

##
# A class containing util methods for database query
class SQLUtil( object ):
    
    # #
    # Search for the Protein in Database having the one of the UniprotAC provided in the list.
    # The protein returned is the one found with he first matching uniprotAC.
    # If several or None protein if found, an exception is raised.
    #
    # @param uniprot_ac_list : string - the list of UniprotAC of the Protein to find
    #
    # @return : Protein - The found Protein
    # @exception : RainetException - If no or several protein is found with the UniprotAC 
    @staticmethod    
    def find_unique_protein_from_uniprotAC_list( uniprot_ac_list ):
        
        if len( uniprot_ac_list ) == 0:
            raise RainetException( "SQLUtil.find_unique_protein_from_uniprotAC_list: Provided uniprotAC list is None." )
        
        for uniprot_ac in uniprot_ac_list:
            protein = SQLUtil.find_unique_protein_from_uniprotAC( uniprot_ac )
            if protein != None:
                return protein
            
        raise RainetException( "SQLUtil.find_unique_protein_from_uniprotAC_list: No Proteins found with uniprotAC list = " + str( uniprot_ac_list ) )
    
    # #
    # Search for the Protein in Database having the provided UniprotAC
    # If several or None protein if found, an exception is raised
    #
    # @param uniprot_ac : string - the UniprotAC of the Protein to find
    #
    # @return : Protein - The found Protein
    # @exception : RainetException - If several protein is found with the UniprotAC 
    @staticmethod
    def find_unique_protein_from_uniprotAC( uniprot_ac ):
        
        sql_session = SQLManager.get_instance().get_session()
        
        protein_list = sql_session.query( Protein ).filter( or_( Protein.uniprotAC == uniprot_ac, Protein.uniprotID == uniprot_ac) ).all()
        # If the protein list is correct (one result), build the relationship between the Protein and the ProteinIsoform
        if protein_list != None and len( protein_list ) != 0:
            if len( protein_list ) == 1:
                protein = protein_list[0]
                if protein != None:
                    return protein
                else:
                    raise RainetException( "SQLUtil.find_unique_protein_from_uniprotAC: Found Proteins with uniprotAC= " + uniprot_ac + " is None." )
            else:
                raise RainetException( "SQLUtil.find_unique_protein_from_uniprotAC: Several Proteins found with uniprotAC= " + uniprot_ac + " : " + str( len( protein_list ) ) + " Proteins" )
        else:
            return None
        
    # #
    # Look at the status of the DB table corresponding to the provided class
    # If the status indicates data were correctly inserted in a previous launch, return True
    # else return False
    #
    # @param table_name : string - The name of the DB table to test
    #
    # @return True if the data in the table were correctly inserted in a previous launch, False if not
    #
    @staticmethod
    def data_already_inserted( table_name ):
        
        # Get a SQL session
        sql_session = SQLManager.get_instance().get_session()
        
        # Look for the status of the given class
        table_status = sql_session.query( TableStatus ).filter( TableStatus.tableName == table_name ).all()
        
        # If there is no entry or several one in the TableStatus table for the provided table, it means the data were not inserted
        if table_status == None or len( table_status ) != 1:
            return False
        
        # Get the status of the entry
        status = table_status[0].tableStatus
        
        # If the table status is OK or WARNING, the data were correctly inserted
        if status == Constants.STATUS_OK or status == Constants.STATUS_WARNING:
            return True
        
        # In other cases, consider the data were not correctly inserted
        return False
    
    
    # #
    # Test whether the given table exists is DB  
    #  
#     @staticmethod
#     def table_exists( table_name ):
#         
#         # Get a SQL session
#         sql_session = SQLManager.get_instance().get_session()
#         
#         # Look for the table in DB
#         try:
#             table_content = sql_session.query( table_name ).all()
#             if table_content == None:
#                 return False
#         except NoSuchTableError:
#             return False
#         
#         return True
    
    # #
    # Clean the table content if required( forced or previous data insertion with error)
    #
    # @param table_name: string - The name of the table to clean
    # @param force_override : boolean - Indicates if the table must be cleaned in all cases
    # @param table_extension : string - (optional) The possible table extension name used in the TableStatus entries
    #
    #
    @staticmethod
    def insert_data_required( table_name, force_override, table_extension = "" ):
              
        if force_override or not SQLUtil.data_already_inserted( table_name + table_extension ):
            return True
        else:
            return False
            
