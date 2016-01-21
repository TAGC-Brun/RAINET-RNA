# ========================================================================================
# This class divides a multiple sequences fasta file in many files as the sequences are
# ========================================================================================



import subprocess as subp
from core.util.format.SeqToFasta import SeqToFasta
from core.util.format.InfoFasta import InfoFasta
from core.util.file.FileUtils import FileUtils
from core.util.log.Logger import Logger



class SplitSeq():
    
    # split_seq
    # ------------
    #
    # split a multiple fasta file in many file
    #
    # Arguments:
    #    - file_sequences_path: file containing the sequences that you want split
    #    - path_ouput: path of output folder
    #    - start header and end header: are the indexes corresponding respectively to first and last identifier letter of header
    #      that will be the name of split file
    # Returns:
    #    - As many fasta file as are the sequence in fasta file
    # 
    # Example:
    # Ensembl header --> '>ENSG00000242389|ENST00000358944|ENSP00000351823'
    # protein id: start = 33, end= 48
    #
    # Uniprot header --> '>sp|P04004|VTNC_HUMAN Vitronectin OS=Homo sapiens GN=VTN PE=1 SV=1'
    # AC: start= 4, end 10
    #
    @staticmethod
    def split_seq(file_sequences_path, path_output, start_header, end_header):
        
        # Through the subprocess method the grep unix command gets the header of fasta file
        # 
        header_dataset = subp.check_output(['grep', '>' , file_sequences_path])
        header = header_dataset.split('\n')
        file_seq = FileUtils.open_text_r(file_sequences_path)
        seq_obj = file_seq.readlines()
            
        for i, term in enumerate(header):
            prot = term[start_header:end_header]
            Logger.get_instance().info(str(i+1) + ' '  + prot)
            
            # extraction of sequence from fasta file
            prot_seq = InfoFasta.get_seq(seq_obj,prot)
            
            # writing of sequence in a fasta file
            fasta_seq = SeqToFasta.give_fasta(term, prot_seq)
            file_out = FileUtils.open_text_w(path_output + prot + '.fasta')
            file_out.write(fasta_seq)
        file_out.close()



