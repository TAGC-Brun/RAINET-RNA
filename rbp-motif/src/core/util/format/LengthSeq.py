#===============================================================================
#===============================================================================
# FIND THE LONGEST SEQEUNCES OF ONE OR MORE GENE FROM FASTA FILE 
#===============================================================================
#===============================================================================


import pickle
from core.util.format.InfoFasta import InfoFasta
from core.util.format.SeqToFasta import SeqToFasta
from core.util.file.FileUtils import FileUtils
from core.util.log.Logger import Logger


class LengthSeq():

    #
    # longest_seq
    # ---------------------------------
    #
    # This method aims to find the longest sequence among isoforms of gene if there are 
    # more sequences with the same length carries out some controls explained below
    # 
    # Arguments:
    #    - seq_obj:file path containing the fasta sequences or list containing the lines of fasta file 
    #    - dict_identifier: file dictionary path that containing the correspondences between the genes and isoforms
    #      this dictionary must be done through the InfoFasta.make_dictionary method giving as input the same seq_obj 
    #      that will be gived as input in this method
    #    - path_ouput_longest: file that will contain the longest seq of gene
    #    - path_output_isoform: file that will contain the isoform with same length
    #    - type_obj: type of seq_obj
    # Returns:
    #    - file with longest seq
    #    - file with isoform seq ----> That will have to be checked (manually or selecting one randomly sequence for gene)
    # 
    # NOTE: 
    # The writing of sequences occurs through SeqToFasta.give_fasta method that returns the sequence in string format
    # this allows to write directly the string in a file
    #
    @staticmethod
    def longest_seq(seq_obj, dict_identifier, path_outfile_longest, path_outfile_isoform, type_obj='list'):
        
        if type_obj == 'file':
            type_text = False
        elif type_obj == 'list':
            type_text = True


        fileout_longest = FileUtils.open_text_a(path_outfile_longest)
        fileout_isoform = FileUtils.open_text_a(path_outfile_isoform)
        file_dict = open(dict_identifier, 'r')
        dict_ids = pickle.load(file_dict)
                               
        
        # Possible conditions:
        #
        # 1) the gene has one longest protein
        #     - in this case this seq is added to longest file 
        # 2) the gene has two protein with the same length
        #    a) the sequences are identical
        #        - in this case one of these identical sequences is added to longest file
        #    b) the sequences are different
        #        - in this case the isoform sequences are added to isoform file
        # 3) the gene has more than two protein with the same length
        #    a) the sequences are identical 
        #        - in this case one of these is added to longest file
        #    b) the sequences are not identical
        #        - the different isoform are added to isoform file
        # 
        # count variables have been initialized in order to check the output during the method elaboration
        #  
        seq_count = 0
        double_seq_count = 0
        not_same_seq_count = 0
        same_seq_count = 0
        more_prot_count = 0
        prot_longest = []
        prot_double_lseq = []
        prot_double_prot = []
        more_two_prot = []
        more_two_lseq = []
        y = 0
        # 
        # This for loop flows on the keys of dictionary 
        for gene in dict_ids:
            y=y+1
            Logger.get_instance().info(str(y) + ' Gene analysed : ' + gene)
            seqs = [] # will contain the isoform list of gene selected
            lenseq = [] # will contain the length of each isoform seq 
            headers = [] # will contain the header of each isoform seq
            
            # this for loop flows on the isoforms of gene selected
            for prot in dict_ids[gene]:
                
                # This lines call InfoFasta class in order to extract
                # the seq, the length and the header of protein selected
                # all item are memorized in lists 
                lenseq.append(InfoFasta.get_length(seq_obj, prot))
                seqs.append(InfoFasta.get_seq(seq_obj, prot, text=type_text))
                headers.append(InfoFasta.get_header(seq_obj,header_identifier=prot, text=type_text))
            
            # Find the max length among the sequences
            # the index_max list contains the index in correspondence of sequence with max length
            len_max = max(lenseq)
            index_max = [item for item in range(len(lenseq)) if lenseq[item] == len_max]
            #
            # The following if conditions check the length of index_max vector
            #
            # Condition 1)
            # -------------
            # if the length of index_max vector is equal to 1 it means that there is just one longest protein
            # the protein sequence is written into longest file 
            #                   
            if len(index_max) == 1:
                Logger.get_instance().info(' If condition 1')
                seq_count +=1
                seq = SeqToFasta.give_fasta(headers[index_max[0]], seqs[index_max[0]]) # (See NOTE above) 
                fileout_longest.write(seq)
                prot_longest.append(dict_ids[gene][index_max[0]])
            #
            # Condition 2)
            # -------------
            # if length of index_max is equal to 2 it means that there are two protein with same length
            #    
            elif len(index_max) == 2:
                Logger.get_instance().info('If condition 2')
                double_seq_count +=1
                # Condition 2a
                # -------------
                # The proteins have the same sequences
                # One protein sequence is written into longest file
                if seqs[index_max[0]] == seqs[index_max[1]]:
                    Logger.get_instance().info('2a')
                    same_seq_count += 1 
                    d_seq = SeqToFasta.give_fasta(headers[index_max[0]], seqs[index_max[0]]) # (See NOTE above)
                    fileout_longest.write(d_seq)
                # Condition 2b
                # -------------
                # The protein have different sequences
                # The sequences are written into isoform file
                else:
                    Logger.get_instance().info('2b')
                    not_same_seq_count += 1
                    for i in range(len(index_max)):
                        prot_double_lseq.append(seqs[index_max[i]])
                        prot_double_prot.append(dict_ids[gene][index_max[i]])
                        diff_seq = SeqToFasta.give_fasta(headers[index_max[i]], seqs[index_max[i]]) # (See NOTE above)
                        fileout_isoform .write(diff_seq)
                        
            # Condition 3)
            # -------------
            # if the length of index_max is greater than two it means that there are more than two proteins 
            # with same length               
            else:
                more_prot_count += 1
                Logger.get_instance().info(' If condition 3')
                
                # Condition 3a
                # -------------
                # The isoforms with same length have actually the same sequences
                # One of this protein is written in longest file
                if seqs.count(seqs[index_max[0]]) == len(index_max):
                    Logger.get_instance().info('3a')
                    m_seq = SeqToFasta.give_fasta(headers[index_max[0]], seqs[index_max[0]]) # (See NOTE above)
                    fileout_longest.write(m_seq)
                    
                # Condition 3b
                # -------------
                # Among the isoforms there are at least two isoforms with different sequences
                # 
                else:
                    Logger.get_instance().info( '3b')
                    more_two_prot.append(gene)
                    more_two_seqs = [] # will contains only the sequences with max length
                    for n in index_max:
                        more_two_seqs.append(seqs[n])
                    more_two_lseq.append(list(set(more_two_seqs))) 
                    for seq in set(more_two_seqs): # set(more_two_seqs) contains only the different sequences
                        # find the sequence index in the list of sequences 
                        index_seq = seqs.index(seq)
                        mdiff_seq = SeqToFasta.give_fasta(headers[index_seq], seqs[index_seq]) # (See NOTE above)
                        fileout_isoform.write(mdiff_seq)
        fileout_isoform.close()
        fileout_longest.close()
        
 
                    
                    


