# This class makes some operations on fasta file

from core.util.file.FileUtils import FileUtils
from core.util.log.Logger import Logger

class InfoFasta():
    
    #
    # get_header
    # -----------
    # 
    # This method get all or one header of a fasta file 
    # Arguments: 
    #    - seq_obj: file path containing the fasta sequences or list containing the lines of fasta file 
    #    - header_identifier: if this arg is specified returns only one header
    #    - text = if the seq_obj is a file path must be False in the other case True
    # Returns: 
    #    - list of header or the header corresponding to header_idenifier
    #
    @staticmethod
    def get_header(seq_obj, header_identifier=None, text=True):
        if text == False:
            fasta = FileUtils.open_text_r(seq_obj)
        elif text == True:
            fasta = seq_obj
        HEADER = []  
        for line in fasta:                      
            if line[0:1] == '>':
                line = line.strip()
                HEADER.append(line)
                

        if header_identifier == None:
            return HEADER
        else:                
            for item in HEADER:
                if header_identifier in item:
                    return item
    
    # 
    # get_seq
    # ------------
    # 
    # This function get the sequence of a seq_obj of a identifier gived as input
    # Arguments:
    #    - seq_obj: file path containing the fasta sequences or list containing the lines of fasta file 
    #    - header_identifier: if this arg is specified returns only one header
    #    - text: if the seq_obj is a file path must be False in the other case True
    # Returns:
    #    - finaleseq: seq selected in string format
    @staticmethod
    def get_seq(seq_obj, header_identifier, text=True):
        if text == False:
            fasta = FileUtils.open_text_r(seq_obj)
        elif text == True:
            fasta = seq_obj
        seq = ''
        flag = 0
        for n, line in enumerate(fasta):
            line = line.strip()
            if line[0:1] == '>' and flag == 0:
                if header_identifier in line:
                    flag = 1
            elif line[0:1] != '>' and flag == 1:
                seq+=line
            elif line[0:1] == '>' and flag == 1:
                finalseq = seq
                flag = 2
            finalseq = seq
        # The replace function is used in order to remove the star because the Ensembl sequences
        # show the star at the end of sequence but if the sequences doesn't show the start anything happens
        return finalseq.replace('*','')
    
    #
    # get_length
    # -------------
    #
    # This function get the length of fasta sequence in a list object   
    #
    # Arguments:
    #    - seq_obj: list containing the lines of fasta file 
    #    - header_identifier: if this arg is specified returns only one header
    # Returns:
    #    - length: length of sequence
    #
    @staticmethod
    def get_length(seq_obj,header_identifier):
        seq = InfoFasta.get_seq(seq_obj, header_identifier, text=True)
        length = len(seq)
        return length
    
    #
    # make_dictionary 
    # ----------------
    # 
    # this method is specific for a fasta file that have the following header format
    # >GENEid|TRANSCRIPTid|PROTEINid
    #
    # This method creates a dictionary with the header of a fasta file
    # The key dictionary is a gene the value is a list of isoform of gene 
    #
    @staticmethod
    def make_dictionary(namefile):
        headers = InfoFasta.get_header(namefile)
        dict_header = {}
        for line in headers:
            line = line[1:]
            ids = line.split('|')
            if len(ids) != 2:
                if ids[0] in dict_header:
                    dict_header[ids[0]].append(ids[2])
                else:
                    dict_header[ids[0]] = [ids[2]]
            else:
                Logger.get_instance().info( " This gene has a protein sequence unavailable : " + line)
                pass
        
        return dict_header
        
        
        
        
        
        
