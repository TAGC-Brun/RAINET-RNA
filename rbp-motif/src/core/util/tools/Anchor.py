# ==============================================================
# This class performs the Anchor analysis calling the local tool
# =============================================================


from core.util import Constants
import subprocess as subp
from core.util.option.OptionManager import OptionManager
from core.util.file.FileUtils import FileUtils
from core.util.log.Logger import Logger
from core.util.option import OptionConstants
from core.util.parsing.FileParser import FileParser
from core.util.parsing.TableWrapper import TableWrapper
from subprocess import CalledProcessError


    


class Anchor():
    
    def __init__(self, output_path):
        
        # Definition of output path
        self.path_output = output_path
        
        # Anchor path definition, it must be a environment variable
        self.anchor_path = Constants.PATH_ANCHOR

    
    # anchor_analysis
    # ----------------------
    # 
    # This method carries out the analysis of Anchor tool (i.e. one only protein)
    # 
    # Arguments:
    #    - fastafile: filename
    #    - motiflist: file containing the motif it will be gived to Anchor tool
    #    - prot: protein identifier
    # Returns:
    #    - Predicted_binding_regions file --> PBR_Prot.txt
    #    - Found_motifs --> Fmotifs_Prot.txt
    #    - Prediction_profile --> Pred_Prot.txt
    #
    def anchor_analysis(self,fastafile, motifslist, prot):
        
        # Calling of anchor command
        # anchor out contains the anchor output in text format
        anchor_out = subp.check_output(["perl", self.anchor_path + 'anchor.pl', fastafile, "-m",  motifslist])
        
        # Definition of the section index of anchor output in order to get a specific section of anchor output
        # Thereby in the next step it will be possible to write these sections in some file
        index_bind_reg = anchor_out.index('Predicted binding regions')
        index_motifs = anchor_out.index('Motifs')
        index_pred_profile = anchor_out.index('Prediction profile output')
        
        # The Anchor output can lack filtered regions section 
        if 'Filtered regions' in anchor_out:
            index_filt_reg = anchor_out.index('Filtered regions')
        else:
            index_filt_reg = index_motifs        

        # ===============================================================================
        # Files writing
        # ===============================================================================
        #
        # Prediction binding regions file ( PBR_protname.txt)
        # This section selects the Prediction binding region of anchor output
        # The PBR section is split in lines and the '#\t' character is removed
        # 
        pbr_out = anchor_out[index_bind_reg:index_filt_reg] 
        pbr_out_list = pbr_out.split('\n')
        pbr = [line[2:] for line in pbr_out_list if line[0:2]=='#\t']
        #
        # When a protein lacks predicting binding region in the output file is written "None" then
        # If the predicting binding regions are not in pbr_out the file writing is skipped 
        #
        if 'None' in pbr_out:
            Logger.get_instance().info( "This protein doesn't contain predicted binding region  ("  + prot +')')
            pass
        elif 'None' not in pbr_out:
            new_pbr_out = [ line.replace(' ', '') for line in pbr]        
            pbr_file = FileUtils.open_text_w(self.path_output +'PBR_' + prot + '.txt')
            pbr_file.write('\n'.join(new_pbr_out))
            pbr_file.close()
        #
        # Found Motifs file (FMotifs_protname.txt)
        #
        fmotifs_out = anchor_out[index_motifs:index_pred_profile]
        fmotifs_out_list = fmotifs_out.split('\n')
        fmotifs = [line[2:] for line in fmotifs_out_list if line[0:2]=='#\t']
        #
        # When a protein lacks Motif in the output file is written "None" then
        # If the Motif are not in fmotif_out the file writing is skipped 
        #
        if 'None' in fmotifs_out:
            Logger.get_instance().info( "This protein doesn't contain any motifs ("  + prot +')')
            pass
        elif 'None' not in pbr_out:
            new_fmotifs = [line.replace(' ','') for line in fmotifs]            
            fmotifs_file = FileUtils.open_text_w(self.path_output+'FMotifs_'+ prot +'.txt')
            fmotifs_file.write('\n'.join(new_fmotifs))
            fmotifs_file.close()
        #
        # Prediction profile output (Pred_protname.txt)
        # This section is always present in anchor output
        #
        pred_file = FileUtils.open_text_w(self.path_output+'Pred_'+ prot+'.txt')
        pred_out = anchor_out[index_pred_profile:]
        string = '#   Columns:\n#   1 - Amino acid number\n#   2 -\
 One letter code\n#   3 - ANCHOR probability value\n#   4 - ANCHOR output\n#'
        pred_out = pred_out.replace(string, '')
        pred_out_list = pred_out.split('\n')
        new_pred_out = [ line.replace(' ', '') for line in pred_out_list]
        final_out = '\n'.join(new_pred_out)
        pred_file.write(final_out)
        pred_file.close()
        
        
    #
    # global_anchor_analysis
    # --------------------------------
    # 
    # This method performs the Anchor analysis for a multiple set of proteins     
    #
    # Arguments:
    #    - file_motifs: file containing the motif it will be given to Anchor tool
    #    - input_folder: directory containing the files to be analyzed
    #    - output_path : directory that will contain all output files
    #
    @staticmethod
    def global_anchor_analysis(file_motifs, input_folder, output_path):
        
        # Description of execution
        Logger.get_instance().info( 'Starting of Anchor Analysis')
        
        # The list file is provided calling a unix command
        try:
            LIST_FILE = subp.check_output(['ls', input_folder])
            LIST_FILE = LIST_FILE.split('\n')
            if '' in LIST_FILE:
                LIST_FILE.remove('')
                for fastafile in LIST_FILE:
                    prot = fastafile.split('.fasta')[0]
                    Logger.get_instance().info(' Anchor Analysis: ' + prot)
                    file_input = input_folder + fastafile
                    # Anchor tool
                    A = Anchor( output_path)
                    A.anchor_analysis( file_input, file_motifs,  prot)
        except CalledProcessError as cpe:
            Logger.get_instance().error( 'Anchor.global_anchor_analysis: Unable to execute listing of files in '+ input_folder)
            raise RbpmotifException( 'Anchor.global_anchor_analysis: Unable to execute listing of files in '+ input_folder, cpe)
        
        Logger.get_instance().info( " End of Anchor Analysis")
        
        
    #
    # anchor_info
    # --------------------
    # 
    # This method reads a Pred_protname.txt file and extracts different information
    # 
    # Arguments:
    #    - input_file: namefile to be analyze
    #    - prot_id: protein name 
    #    - nume_aa: number of amino acids that will be taken in account in order to find the up regions
    # Returns:
    #    - dictionary: containing anchor informations
    #
    @staticmethod
    def anchor_info(input_file, prot_id, num_aa=6):
        
        Logger.get_instance().info( ' Creation of a dictionary containing the dictionary analysis \
informations of protein ' + prot_id ) 
        
        
        # Initialization of a dictionary that will contain the output informations 
        dictionary_anchor = {}
        
        # Reading of input file containing the anchor output for one protein 
        anchor_table = FileParser.make_table(input_file, '\t', skip=2)
        
        # Extraction of anchor output information
        position = TableWrapper.get_column(anchor_table, 0)
        aminoacid = TableWrapper.get_column(anchor_table, 1)
        probability = TableWrapper.get_column(anchor_table, 2)
        binary_array = TableWrapper.get_column(anchor_table, 3)
        
        # Conversion in float number the probability array
        prob = [ float(item) for item in probability]
        # Conversion in int number of binary array 
        binary_array = [ int(item) for item in binary_array]
        
        # Storing informations
        dictionary_anchor["position"] =  position
        dictionary_anchor["aminoacids"] = aminoacid
        dictionary_anchor["probability"] = prob
        
        # Counting of ones number in binary array 
        # This represent the number of aminoacids that has a probability prediction greater than 0.5
        num_ones = binary_array.count(1)
        length = len(binary_array)
        # Fraction calculation
        fraction = round(num_ones/float(length), 3)
        
        dictionary_anchor["values>threshold"] = num_ones
        dictionary_anchor["fraction"] = str(fraction)
        #
        # This section counts the regions that have almost 6 consecutive aminoacid (Anchor --> min length 6 AA)
        # for each region found the start and end positions are taken and memorized in a vector up_region
        count_one = 0
        up_region = []
        #
        # This for loop flows over the binary array and checks if the i-th aminoacid number is one or zero
        # if it is a one the variable count_one is increased by one
        #
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
                if count_one >= num_aa:
                    start = end - count_one
                    region = [start, end]
                    count_one = 0
                    up_region.append(region)
                #
                # if the count one is not >= 10 the count_one is reseted
                else:
                    count_one = 0
            #
            # this section is necessary when the last number of binary array is not a zero
            # indeed it would not be possible to memorized the last ones region
            #
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
        
        # The informations are stored into dictionary 
        
        dictionary_anchor["binary_array"] = binary_array
        dictionary_anchor["regions"] = up_region
        dictionary_anchor["num region"] = str(num_region_up)
        dictionary_anchor["anchor table"] = table

        return dictionary_anchor


    # anchor_string_info
    # -------------------------------
    #
    # This method returns three different string in order to create a file
    # 
    # Arguments: 
    #    - input_file: file name that you want analyze
    #    - prod_id: protein name 
    #    - num_aa: number of amino acids that will be taken in account in order to find the up regions
    #
    # Returns:
    #     string 1 : protein_name, fraction, num_region
    #     string 2 : protein_name, N, Start, End, Length_region 
    #     (the second file is provided directly also by ANCHOR tools but in order to have a control this method returns it all the same)
    #
    @staticmethod
    def anchor_string_info(input_file, prot_id, num_aa):
        
        # This line call the previously method in order to find the ones region
        # all informations are contained in a dictionary
        out = Anchor.anchor_info(input_file, prot_id, num_aa=num_aa)
        
        
        fraction = out["fraction"]
        num_region = out["num region"]

        # The following lines allow to create the strings containing some
        # dictionary informations that will can be put into a file
        
        row_file_1 = [prot_id, fraction, num_region]
        string_file_1 = '\t'.join(row_file_1)        
        
        
        row_file = out["anchor table"]
        if row_file != []:
            lines_file_2 = ['\t'.join(line) for line in row_file]        
            string_file_2 = '\n'.join(lines_file_2)
            
        # When the row_file is empty it means that this protein haven't any up_region
        # Then in order to not lose the protein information a row of zeros is added
        elif row_file == []:
            row = [prot_id, '0', '0', '0', '0']
            string_file_2 =  '\t'.join(row)

        # With this method are created two text type (see the lines before method definition)
        string_file = ( string_file_1, string_file_2)
        
        return string_file
    
    # make_anchor_file
    # --------------------
    #    
    # This method carries out the previously methods for one or more protein 
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
    def make_anchor_file(input_path, output_path, num_aa, dataset_type):
        
        
        # initialization of file names
        file_name_1 = dataset_type + '_AnchorTable' + '.txt'
        file_name_2 = dataset_type + '_AnchorRegion' + '.txt'

        
        num_aa_string = '('+ str(num_aa) +' AA)'
        
        # Files opening and title string writing
        file_1 = FileUtils.open_text_a(output_path + file_name_1)
        file_2 = FileUtils.open_text_a(output_path + file_name_2)
       
        header_file_table = ['Protein', 'Fraction', 'Region N.'+ num_aa_string]
        header_file_region = ['Protein', 'N', 'Start', 'End',  'Region length']
        
        header_string_table = '\t'.join(header_file_table)
        header_file_region = '\t'.join(header_file_region)
        
        file_1.write(header_string_table + '\n')
        file_2.write(header_file_region + '\n')

        # This command allows to taken the file names of protein that you want analyze
        list_file = subp.check_output(['ls', input_path])
        list_file = list_file.split('\n')
        if '' in list_file:
            list_file.remove('')


        # This section performs the iupred_string_info method (that calls also iupred_info method) 
        # for each protein file in list_file and simultaneously appends into files the output results
        # in a tab format
        i = 0
        for n, pred_file in enumerate(list_file):
            if pred_file.split('_')[0] == 'Pred':
                i += 1         
                prot_id = pred_file.split('.')[0].split('_')[1]
                Logger.get_instance().info( str(i) + ' ' + prot_id)
                namefile = input_path + pred_file
                out_string = Anchor.anchor_string_info(namefile, prot_id, num_aa)
                
                string_file_1 = out_string[0]
                string_file_2 = out_string[1]
                
                file_1.write(string_file_1 + '\n')
                file_2.write(string_file_2 + '\n')

            
            
        file_1.close()
        file_2.close()
            
  
if __name__ == '__main__':
    
    OptionManager.get_instance().initialize()
    
    Anchor.global_anchor_analysis(OptionManager.get_instance().get_option( OptionConstants.OPTION_FILENAME),
                                OptionManager.get_instance().get_option( OptionConstants.OPTION_INPUT_PATH),
                                OptionManager.get_instance().get_option( OptionConstants.OPTION_OUTPUT_PATH))

