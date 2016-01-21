# =============================================================
# This class allows to access to Uniprot through URL connection
# =============================================================


from core.util.web.Web import Web
from core.util.parsing.FileParser import FileParser
from core.util.file.FileUtils import FileUtils
from core.util.log.Logger import Logger
from core.util import Constants


class Uniprot():
    
    def __init__(self):
        
        self.server = Constants.SERVER_UNIPROT
        self.end = Constants.END_URL_UNIPROT
    #   
    # connection
    # ----------------
    #   
    # This method allows the connection to Uniprot database
    #
    # Arguments:
    #    -qeury: can be a protein id or an accassion number
    #    -type_query: can be fasta or info
    # Returns:
    #    -text page containing the sequence or info protein
    #
    def connection(self, query, type_query): 
        
        if type_query == 'fasta':
            end = self.end[1]
        elif type_query != 'fasta':
            end = self.end[0]
        url = self.server + query + end
        text_page = Web.openurl(url)
        
        return text_page
    
    #
    # get_sequence
    # ----------------
    #
    # This method gets the sequences of query protein
    # format_out = True --> string fasta format
    # format_out != True --> sequence in string 
    #
    @staticmethod 
    def get_sequence(query, format_out=True): 
        
        type_query = 'fasta'
        U = Uniprot()
        fasta_seq = U.connection(query, type_query)
        
        if format_out == True:
            return fasta_seq
        else:
            string_seq = []
            for line in fasta_seq:
                line = line.strip()
                string_seq.append(line)
    #
    # get_info
    # ----------------
    #
    # This method gets the info of query protein
    #            
    @staticmethod
    def get_info(query):
        
        type_query = 'info'
        U = Uniprot()
        info_protein = U.connection(query, type_query)
      
        
        return info_protein
    
    #
    # get_list_seq
    # ------------------
    # 
    # This method gets the sequences of a protein list (path_input_list)
    # returns the sequences in a file (path_output)
    #  
    @staticmethod
    def get_list_seq(path_input_list, path_output):
        
        seq_file = FileUtils.open_text_a(path_output)
        
        
        protein_list = FileParser.read_file(path_input_list)
        for protein in protein_list:
            seq = Uniprot.get_sequence(protein, format_out=True)
            seq_file.write(seq )
            
        seq_file.close()
            
            
        
        
        
        
                
        
        
        