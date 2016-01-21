# =========================================================================
# This class performs the RNA-binding analysis of DisoRDPbind output
# =========================================================================


from core.util import Constants
from core.util.file.FileUtils import FileUtils
from core.util.log.Logger import Logger



class DisoRDPbind():
    
    # output_reading
    # -------------------
    # 
    # reading of DisoRDPbind output file
    # Argument:
    #    - filename
    # Returns:
    #    - list of string: each item containing one protein output 

    @staticmethod
    def output_reading(filename):
        
        input_file = FileUtils.open_text_r(filename)
        
        text_file = []
        
        lines = input_file.readlines()
        string = ''
        for n, line in enumerate(lines):
            if line[0:1] == '>' and n == 0:
                string += line[1:] 
            elif line[0:1] != '>' and n!=0:
                string += line 
            elif line[0:1] == '>' and n!=0:
                # append in string format the output of one protein
                text_file.append(string)
                # reset the string variable and add the header
                string = ''
                string += line[1:]
            else:
                Logger.get_instance().info( ' Check this line : ' +line)


        text_file.append(string)        
        
        return text_file
    
    
    # fraction_calculation
    # --------------------------
    # 
    # this method counts the number of residues that are predicted as putative disordered RNA(-DNA,-protein) binding residues
    # this method takes as input the output of one protein in string format
    # The sum of putative disordered X-binding residues is normalized to length of sequences
    #
    # Arguments:
    #    - protein_output: string of one output protein
    #    - binding partner: identify the prediction type that you want analyzed 
    #         1----> RNA-binding
    #         2----> DNA-binding
    #         3----> protein-binding
    #    - num_aa: number of amino acids that will be taken in account in order to find the up regions
    #
    # Returns:
    #    - dictionary containing the output informations
    #
    @staticmethod
    def fraction_calculation(protein_output, binding_partner, num_aa=10):
        
        # the protein output is split in lines
        protein_lines = protein_output.split('\n')
        
        # output format
        # line0 ---> protein name
        # line1 ---> sequence
        # line2 ---> RNA-binding residues
        # line3 ---> RNA-binding prediction
        # line4 ---> DNA-binding residues
        # line5 ---> DNA-binding prediction
        # line6 ---> protein-binding residues
        # line7 ---> protein-binding prediction
        
        
        protname = protein_lines[0]
        sequence = protein_lines[1]
        aminoacid = [letter.upper() for letter in sequence]
        length_seq = len(sequence)
        
        Logger.get_instance().info( ' Creation of a dictionary containing the DisoRDPbind analysis \
informations of protein ' + protname )
        
        # this vector represent the position of sequence
        position_array = range(1,length_seq + 1)
        position  = [str(item) for item in position_array]
        #
        # this dictionary contains the keyword in order to extract a specific type of binding residues        
        #
        dict_bind_partner = Constants.BINDING_PARTNER
        #
        # this line gets the line index corresponding to the desired binding residues 
        #
        index_binding_partner = [protein_lines.index(line) for line in protein_lines if dict_bind_partner[binding_partner] in line]
        binary_array_char = protein_lines[index_binding_partner[0]].split(':')[1]
        probability = protein_lines[index_binding_partner[0] + 1].split(':')[1]
        #
        # this is the binary array containing the information of putative disordered residues 
        binary_array = [int(item) for item in binary_array_char]
        #
        # number of putative disordered residues
        #
        num_ones = binary_array.count(1)
        length = len(binary_array)
        fraction = round(num_ones/float(length), 3)
        #
        # Information Storing 
        #
        dictionary_rdpbind = dict([('protein name', protname), ('position', position), ('sequence', aminoacid), ('probability', probability), ('values>threshold', num_ones),
                                  ('fraction', str(fraction)), ('binary array', binary_array)])
        
        
        
        #
        # This section counts the regions that have almost 10 consecutive aminoacid
        # for each region found the start and end positions are taken and memorized in a vector up_region
        #
        count_one = 0
        up_region = []
        #
        # This for loop flows over the binary array and checks if the i-th aminoacid number is one or zero
        # if it is a one the variable count_one is increased by one 
        #
        for ind, num in enumerate(binary_array):
            if num == 1:
                count_one+=1
            # when the for loop bumps into a zero the ones counting is stopped and  
            # the count_one is compared with number of 10 that representing the 10 aminoacids            
            elif num == 0:
                end = ind 
                # if the count_one is effectively >= 10 the start and end positions of region are 
                # memorized in the up_region vector
                if count_one >= num_aa:
                    start = end - count_one
                    region = [start, end]
                    count_one = 0
                    up_region.append(region)
                # if the count one is not >= 10 the count_one is reseted
                else:
                    count_one = 0
            #
            # this section is necessary when the last number of binary array is not a zero
            # indeed it would not be possible to memorized the last ones region
            # this if condition will be true when the binary array shows an one in the last i-th number
            # in this way it is possible memorized the last ones region in the vector up_region
            #
            if num == 1 and ind == len(binary_array) - 1:
                end =  ind + 1
                if count_one >= num_aa:
                    start = end - count_one
                    region = [start, end]
                    count_one = 0
                    up_region.append(region)
                else:
                    count_one = 0
            else:
                pass
        #
        # This part creates a table with specific positions of ones regions
        # indeed the start and end position of up_region vector are mapped over the position vector
        # in order to extract the exactly positions
        #
        table = []
        for n, region in enumerate(up_region):
            row = []
            row.append(protname)
            n +=1
            numb = str(n)
            row.append(numb)
            start = position[region[0]]
            end = position[region[1] - 1]
            row.append(start)
            row.append(end)
            length = region[1] - region[0]
            row.append(str(length))
            table.append(row)
                    
            
        num_region_up = len(up_region)
        
        
        dictionary_rdpbind["regions"] = up_region
        dictionary_rdpbind["num region"] = str(num_region_up)
        dictionary_rdpbind["DisoRDPbind table"] = table
        
        return dictionary_rdpbind
        
    
    # disordp_string_info
    # ----------------------------
    # 
    # This method returns two different string in order to create a file
    # 
    # Arguments: 
    #    - protein_output: string of output protein that you want analyze
    #    - binding_partner: see previously method
    #    - num_aa: number of amino acids that will be taken in account in order to find the up region
    #
    # Returns:
    #     string 1 : protein_name, fraction, num_region
    #     string 2 : pritein_name, N, Start, End, Length_region
    #
    @staticmethod
    def disordp_string_info(protein_output, binding_partner, num_aa):
        
        
        # This line call the previously method in order to find the ones region
        
        out_disordp = DisoRDPbind.fraction_calculation(protein_output, binding_partner, num_aa)
        
        
        # Information extraction
        
        protname = out_disordp["protein name"]
        fraction = out_disordp["fraction"]
        num_region = out_disordp["num region"]

        Logger.get_instance().info("Protein analysed: "+ protname)
        
        # The following lines allow to create the strings containing some
        # dictionary informations that will can be put into a file
        
        row_file_1 = [protname, fraction, num_region]
        string_file_1 = '\t'.join(row_file_1)
        
        
        
        row_file_2 = out_disordp["DisoRDPbind table"]
        if row_file_2 != []:
            lines_file_2 = ['\t'.join(line) for line in row_file_2]        
            string_file_2 = '\n'.join(lines_file_2)
        # When the row_file is empty it means that this protein haven't any up_region
        # Then in order to not lose the protein information a row of zeros is added
        #
        elif row_file_2 == []:
            row_2 = [protname, '0', '0', '0', '0']
            string_file_2 =  '\t'.join(row_2)
           
        

        # With this method are created two text type (see the lines before method definition)
        string_file = ( string_file_1, string_file_2)
        
        
        return string_file


    # make_disorbd_file
    # -----------------------------
    # 
    # This method carries out the previously methods for the protein contained in a DisoRDPbind
    # output file and returns directly two file containing the DisoRDPbind string information 
    #
    # Arguments: 
    #    - input_path: directory path containing the files that you want analyzed
    #    - output_path: output folder
    #    - num_aa: the same of previously methods
    #    - dataset_type: prefix to be add to output file
    #
    # Returns:
    #    - two file containing the anchor string information 
    #

    @staticmethod
    def make_disorbd_file(input_path, output_path, binding_partner, num_aa, dataset_type):
        
        
        # initialization of file names
        file_name_1 = dataset_type + '_DisoRDPbindTable.txt'
        file_name_2 = dataset_type + '_DisoRDPbindRegion.txt'
        
        
        num_aa_string = '('+ str(num_aa) +' AA)'
        
        # Files opening and title string writing
        file_1 = FileUtils.open_text_a(output_path + file_name_1)
        file_2 = FileUtils.open_text_a(output_path + file_name_2)
        header_file_table = ['Protein', 'Fraction ','Region N.' +  num_aa_string]
        header_file_region = ['Protein', 'N', 'Start', 'End',  'Region length']
        
        
        header_string_table = '\t'.join(header_file_table)
        header_file_region = '\t'.join(header_file_region)
        
        file_1.write(header_string_table + '\n')
        file_2.write(header_file_region + '\n')
        
        
        # Reading of DisoRDPbind output file 
        protein_output_list = DisoRDPbind.output_reading(input_path)
        
        for n, output in enumerate(protein_output_list):
            if 'WARNING:' in output:
                prot = output.split('\n')[0]
                Logger.get_instance().info( str(n+1) + "\n This protein contains >=10000 residues\
 (DisoRBDbind cannot predict the proteins with size >=10000) " + prot)
            else:
                Logger.get_instance().info( str(n+1))
                results = DisoRDPbind.disordp_string_info(output, binding_partner, num_aa)
            
        
                string_file_1 = results[0]
                string_file_2 = results[1]
        
                file_1.write(string_file_1 + '\n')
                file_2.write(string_file_2 + '\n')
        
        
        file_1.close()
        file_2.close()
        
        
        
        
        
        
        
    
    
        
    
    