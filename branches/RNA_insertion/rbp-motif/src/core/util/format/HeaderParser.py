    
from core.util.file.FileUtils import FileUtils
from core.util.log.Logger import Logger
from core.util.parsing.FileParser import FileParser
    


class HeaderParser():
    
    
    # change_header
    # ---------------------
    # This function change the header of fasta file
    # 
    # Arguments: 
    #    - path_input_file: fasta fale of which you want change the header
    #    - path_output_file: path of output file
    #    - source: represent the resource from which come the input file
    #         source 1---> Ensembl
    #         source 2---> Uniprot
    #    - type_id: only when source = 2 you can set this parameter:
    #            type_id = 1 ---> AC
    #            type_id = 2---> identifier
    #
    # Note: if the source is Ensembl the new header is composed by protein ID
    #    
    @staticmethod
    def change_header(path_input_file, path_ouptut_file, source=1, type_id=1):
        
        file_input = FileUtils.open_text_r(path_input_file)
        seq_file = file_input.read()

        
        file_output = FileUtils.open_text_a(path_ouptut_file)
        
        # Warning: Check that the file have the '>' char only at beginning of header lines and not in other points
        # otherwise the split will occur in an incorrect way!
        seq_file_list = seq_file.split('>')[1:] 

        
        for seq in seq_file_list:
            lines = seq.split('\n')
            header = lines[0]
            Logger.get_instance().info(header)
            # Ensembl
            if source == 1:
                new_header = '>' + header.split('|')[2] +'\n' # see Note
            # Uniprot
            elif source == 2:
                diff_header = header.split(' ')[0]
                # AC
                if type_id == 1:                    
                    new_header =  '>' + diff_header.split('|')[1] + '\n'
                # ID
                elif type_id == 2:
                    new_header =  '>' + diff_header.split('|')[2] + '\n'
                    
            fasta = new_header + '\n'.join(lines[1:])
            
            file_output.write(fasta)
            
        file_output.close()
             
               
                    
            