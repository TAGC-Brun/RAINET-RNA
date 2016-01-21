# This class generate a string in fasta format starting from one header and sequence(string format)


from math import fmod

class SeqToFasta():
    
    # give_fasta
    # ----------
    #
    # Arguments:
    #    - header: will be the header of sequence
    #    - sequence: must be in string format 
    # Returns:
    #    - fasta : a string containing the fasta format of sequence (60 AA for line)
    #
    
    
    @staticmethod
    def give_fasta( header, sequence):
        fasta ='' 
        fasta += header +'\n'
        # count the number of AA
        length = len(sequence)
        # if the sequence has more than 60 AA
        # count the number of line composed by 60 AA
        if length >= 60:
            row = length/60
            start = 0
            end = 60
            for i in range(row):
                fasta += sequence[start:end] + '\n'
                start = end # the start becomes the end
                end += 60 # the end is increased of 60
            if fmod(length,60) != 0: # if the remainder is not equal to 0 means that there are others AA to be add
                fasta += sequence[start:] + '\n'
        # if the sequence has less than 60 AA
        # just one line will be written
        elif length < 60:
            fasta += sequence + '\n'
        return fasta

        
        
        
        
        


