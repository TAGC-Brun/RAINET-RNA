# ==============================================================
# This class performs the Iupred analysis calling the local tool
# =============================================================


import subprocess as subp
from core.util.option.OptionManager import OptionManager
from core.util.file.FileUtils import FileUtils
from core.util.log.Logger import Logger
from core.util.option import OptionConstants
from core.util import Constants
from core.util.parsing.FileParser import FileParser
from core.util.parsing.TableWrapper import TableWrapper
from subprocess import CalledProcessError
from core.util.exception.RbpmotifException import RbpmotifException



class Iupred():
    def __init__(self, output_path):
        
        # Definition of output path
        self.path_output = output_path 
        
        # IUPred_PATH definition, it must be a environment variable
        self.iupred_path = Constants.PATH_IUPRED

        
    # iupred_analysis
    # ----------------------
    # 
    # This method carries out an analysis of IUPred tool (i.e. one only protein)  
    # 
    # Arguments:
    #    - fastafile: filename
    #    - prot: protein identifier
    # Returns:
    #    - Prediction_profile --> Iupred_Prot.txt
  
    def iupred_analysis(self, fastafile, prot):
        
        self.fastafile = fastafile
        self.prot = prot
        
        # Calling of IUPred command
        iupred_out = subp.check_output([self.iupred_path+"iupred", self.fastafile, "long"])
        
        
        # Prediction output file ( Prediction_protname.txt)
        pred_file = FileUtils.open_text_w(self.path_output+'IUPred_' + self.prot + '.txt')
        index_prediction = iupred_out.index('Prediction output')
        iupred_out = iupred_out[index_prediction:] 
        iupred_out_list = iupred_out.split('\n')
        new_iupred_out = []
        new_iupred_out.append(iupred_out_list[0])
        new_iupred_out.append(iupred_out_list[1])
        for line in iupred_out_list[2:]:
            new_line = []
            for item in line.split(' '):
                if item != '':
                    new_line.append(item)
            new_line_string = '\t'.join(new_line)
            new_iupred_out.append(new_line_string)
                
        final_out = '\n'.join(new_iupred_out) 
            
        pred_file.write(final_out)
        pred_file.close()
        
    #
    # global_iupred_analysis
    # --------------------------------
    # 
    # This method performs the Iupred analysis for a multiple set of proteins     
    #
    # Arguments:

    #    - input_folder: directory containing the files to be analyzed
    #    - output_path : directory that will contain all output files
    #   
    @staticmethod
    def global_iupred_analysis(input_folder, output_path):
        
        # Description of execution
        Logger.get_instance().info( 'Starting of Iupred Analyisis')
        
        try:
            list_file = subp.check_output(['ls', input_folder])
            list_file = list_file.split('\n')
            if '' in list_file:
                list_file.remove('')
            for fastafile in list_file:
                prot = fastafile.split('.fasta')[0]
                Logger.get_instance().info('Iupred Analysis : ' + prot)
                file_input = input_folder + fastafile
                # IUPred tool
                I = Iupred(output_path)
                I.iupred_analysis(file_input, prot)
        except CalledProcessError as cpe:
            Logger.get_instance().error( 'Iupred.global_anchor_analysis: Unable to execute listing of files in '+ input_folder)
            raise RbpmotifException( 'Iupred.global_anchor_analysis: Unable to execute listing of files in '+ input_folder, cpe)
        
                    
                
        Logger.get_instance().info( " End of Iupred Analysis")
    
    #
    # iupred_info
    # --------------------
    # 
    # This method reads a IUPred_protname.txt file and extract different information
    # 
    # Arguments:
    #    - input_file: namefile to be analyze
    #    - prot_id: protein name 
    #    - nume_aa: number of amino acids that will be taken in account in order to find the up regions
    #    - threshold: values of probability considered to be a significant prediction value
    # Returns:
    #    - dictionary: containing anchor informations
    #
    #
    # Keys : position, aminoacid, probability, fraction, values>threshold, regions, binary_array, num_region, iupred_table
    #
    @staticmethod
    def iupred_info(input_file, prot_id, threshold, num_aa=10):
        
        Logger.get_instance().debug( ' Creation of a dictionary containing the iupred analysis \
informations of protein ' + prot_id + '. The threshold adopted is: ' + str(threshold))
        
        # Initialization of a dictionary that will contain the output informations
        dictionary_iupred = {}
        
        # Reading of input file containing the anchor output for on protein 
        iupred_table = FileParser.make_table(input_file, '\t', skip=2)
        
        # Extraction of anchor output information
        position = TableWrapper.get_column(iupred_table, 0)
        aminoacid = TableWrapper.get_column(iupred_table, 1)
        probability = TableWrapper.get_column(iupred_table, 2)
        
        prob = [ float(item) for item in probability]
        
        dictionary_iupred["position"] =  position
        dictionary_iupred["aminoacid"] = aminoacid
        dictionary_iupred["probability"] = prob
        
        
        # this for loop flow over the probability array in order to make a binary array: 
        # if the value is >= threshold set to one number and set a zero if it isn't
        binary_array = []
        for item in prob:
            if item >= threshold:
                binary_array.append(1)
            else:
                binary_array.append(0)
        
        
        num_ones = binary_array.count(1)
        length = len(binary_array)
        fraction = round(num_ones/float(length), 3)
        
        dictionary_iupred["values>threshold"] = num_ones
        dictionary_iupred["fraction"] = str(fraction)
        
        # This section counts the regions that have almost 10 consecutive aminoacid
        # for each region found the start and end positions are taken and memorized in a vector up_region  
        count_one = 0
        up_region = []
        #
        # This for loop flows over the binary array and checks if the i-th aminoacid number is one or zero
        # if it is a one the variable count_one is increased by one 
        for ind, num in enumerate(binary_array):
            if num == 1:
                count_one+=1
            #
            # when the for loop bumps into a zero the ones counting is stopped and  
            # the count_one is compared with number of 10 that representing the 10 aminoacids
            #
            elif num == 0:
                end = ind 
                # if the count_one is effectively >= 10 the start and end positions of region are 
                # memorized in the up_region vector
                #
                if count_one >= num_aa:
                    start = end - count_one
                    region = [start, end]
                    count_one = 0
                    up_region.append(region)
                #
                # if the count one is not >= 10 the count_one is reseted
                #
                else:
                    count_one = 0
            #
            # this section is necessary when the last number of binary array is not a zero
            # indeed it would not be possible to memorized the last ones region
            # this if condition will be true when the binary array shows an one in the last i-th number
            # in this way it is possible memorized the last ones region in the vector up_region
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
        table = []
        for n, region in enumerate(up_region):
            row = []
            row.append(prot_id)
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
        
        #  informations storing
        #
        dictionary_iupred["binary_array"] = binary_array
        dictionary_iupred["regions"] = up_region
        dictionary_iupred["num region"] = str(num_region_up)
        dictionary_iupred["iupred table"] = table
            
        return dictionary_iupred
    
    
    # iupred_string_info
    # -------------------------------
    #
    # This method returns three different string in order to create some files
    # 
    # Arguments: 
    #    - input_file: file name that you want analyze
    #    - prod_id: protein name 
    #    - num_aa: number of amino acids that will be taken in account in order to find the up regions
    #    - threshold 1 and 2 : this method allows to give two values of threshold that will be compared 
    #
    # Returns:
    #     string 1 : protein_name, fraction th_1, fraction th_2, num_region_th1, num_region_th2
    #     string 2 : pritein_name, N, Start, End, Length_region (th1)
    #     string 3 : pritein_name, N, Start, End, Length_region (th2)
    #
    # This method returns three different string in order to create a file
    #
    @staticmethod
    def iupred_string_info(input_file, prot_id, threshold_1, threshold_2, num_aa):
        
        # This line call the previously method in order to find the ones region
        # all informations are contained in a dictionary
        # the method is called twice in order to extract the information for two threshold values
        
        out_th1 = Iupred.iupred_info(input_file, prot_id, threshold_1, num_aa=num_aa)
        out_th2 = Iupred.iupred_info(input_file, prot_id, threshold_2, num_aa=num_aa)
        
        
        fraction_th1 = out_th1["fraction"]
        fraction_th2 = out_th2["fraction"]
        num_region_th1 = out_th1["num region"]
        num_region_th2 = out_th2["num region"]
        
        
        # The following lines allow to create the strings containing some
        # dictionary informations that will can be put into a file
        
        row_file_1 = [prot_id, fraction_th1, fraction_th2, num_region_th1, num_region_th2 ]
        string_file_1 = '\t'.join(row_file_1)
        
        
        
        row_file_2 = out_th1["iupred table"]
        if row_file_2 != []:
            lines_file_2 = ['\t'.join(line) for line in row_file_2]        
            string_file_2 = '\n'.join(lines_file_2)
        # When the row_file is empty it means that this protein haven't any up_region
        # Then in order to not lose the protein information a row of zeros is added
        elif row_file_2 == []:
            row_2 = [prot_id, '0', '0', '0', '0']
            string_file_2 =  '\t'.join(row_2)
           
        
        
        row_file_3 = out_th2["iupred table"]
        if row_file_3 != []:
            lines_file_3 = ['\t'.join(line) for line in row_file_3]        
            string_file_3 = '\n'.join(lines_file_3)
        elif row_file_3 == []:
            row_3 = [prot_id, '0', '0', '0', '0']
            string_file_3=  '\t'.join(row_3)
        

        # With this method are created three text type (see the lines before method definition)
        string_file = ( string_file_1, string_file_2, string_file_3)
        
        
        return string_file
    
    # make_iupred_file
    # --------------------
    #    
    # This method carries out the previously methods for one or more protein 
    # 
    # Arguments: 
    #    - input_path: directory path containing the files that you want analyzed
    #    - output_path: output folder
    #    - num_aa: the same of previously methods
    #    - dataset_type: prefix to be add to output file
    #    - th_1, th_2: threshold values
    #
    # Returns:
    #    - two file containing the anchor string information 
    #

    
    @staticmethod
    def make_iupred_file(input_path, output_path, th_1, th_2, num_aa, dataset_type):
        
        
        # initialization of file names
        file_name_1 = dataset_type + '_IupredTable' + '_t1_'+ str(th_1) + '_t2_' + str(th_2) + '.txt'
        file_name_2 = dataset_type + '_IupredRegion_' + str(th_1)  + '.txt'
        file_name_3 = dataset_type + '_IupredRegion_' + str(th_2)  + '.txt'
        
        
        num_aa_string = '('+ str(num_aa) +' AA)'
        
        # Files opening and title string writing
        file_1 = FileUtils.open_text_a(output_path + file_name_1)
        file_2 = FileUtils.open_text_a(output_path + file_name_2)
        file_3 = FileUtils.open_text_a(output_path + file_name_3)
        header_file_table = ['Protein', 'Fraction '+ str(th_1), 'Fraction ' + str(th_2), 'Region N.' +  num_aa_string +'th_'+ str(th_1) , 'Region N.'+ num_aa_string +'th_'+ str(th_2)]
        header_file_region = ['Protein', 'N', 'Start', 'End',  'Region length']
        
        header_string_table = '\t'.join(header_file_table)
        header_file_region = '\t'.join(header_file_region)
        
        file_1.write(header_string_table + '\n')
        file_2.write(header_file_region + '\n')
        file_3.write(header_file_region + '\n')
        
        # This command allows to taken the file names of protein that you want analyze
        list_file = subp.check_output(['ls', input_path])
        list_file = list_file.split('\n')
        if '' in list_file:
            list_file.remove('')

        # This section performs the iupred_string_info method (that calls also iupred_info method) 
        # for each protein file in list_file and simultaneously appends into files the output results
        # in a tab format
        for i, pred_file in enumerate(list_file):
            i += 1
            prot_id = pred_file.split('.')[0].split('_')[1]
            Logger.get_instance().info( str(i) + ' ' + prot_id)
            namefile = input_path + pred_file
            out_string = Iupred.iupred_string_info(namefile, prot_id, th_1, th_2, num_aa)
            
            string_file_1 = out_string[0]
            string_file_2 = out_string[1]
            string_file_3 = out_string[2]
            
            file_1.write(string_file_1 + '\n')
            file_2.write(string_file_2 + '\n')
            file_3.write(string_file_3 + '\n')
            
            
        file_1.close()
        file_2.close()
        file_3.close()
            
            
        
               
        
        
        
        
        
if __name__ == '__main__':
    
    OptionManager.get_instance().initialize()
    
    Iupred.global_iupred_analysis(OptionManager.get_instance().get_option( OptionConstants.OPTION_INPUT_PATH),
                                OptionManager.get_instance().get_option( OptionConstants.OPTION_OUTPUT_PATH))


        
    