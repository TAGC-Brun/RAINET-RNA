# ===========================================================================================
# This class performs a Disorder Analysis of a dataset with Anchor, Iupred, DisoRDPbind tools
# ===========================================================================================


# Starting files and conditions
#    - fasta file containing the proteins sequences of a certain dataset
#    - Anchor and Iupred tools must be installed in local
#    - motifs list to provide to Anchor tool
#    - the DisoRDPbind output obtained previously by online tool
#



from core.util import Constants
from core.util.property.PropertyManager import PropertyManager
from core.data import DataConstants
from core.util.option import OptionConstants
from core.util.option.OptionManager import OptionManager
from core.util.time.Timer import Timer

from core.util.tools.Anchor import Anchor
from core.util.tools.Iupred import Iupred
from core.util.format.HeaderParser import HeaderParser
from core.util.log.Logger import Logger
from core.util.format.SplitSeq import SplitSeq
from core.util.parsing.FileParser import FileParser
from core.util.parsing.TableWrapper import TableWrapper
from core.util.file.FileWriter import FileWriter
from core.util.tools.DisoRDPbind import DisoRDPbind


class DisorderAnalysis():
    


    # Header change
    # ----------------------------
    #
    # This method modifies the header of a fasta file (if you have need to modify it)
    # This is not a mandatory step
    # 
    # Output: new file in the path_file_output
    #
    def change_header(self):
        
        # Variables definition
        
        self.path_home = Constants.PATH_HOME
        
        self.path_input =  self.path_home +  PropertyManager.get_instance().get_property( DataConstants.HEADER_INPUT_SEQ_PROPERTY, True)
        self.namefile = PropertyManager.get_instance().get_property( DataConstants.HEADER_FILE_SEQ_PROPERTY, True)
        self.path_file_input = self.path_input + self.namefile
        
        self.path_output =  self.path_home +  PropertyManager.get_instance().get_property( DataConstants.HEADER_OUTPUT_SEQ_PROPERTY, True)
        self.path_file_output = self.path_output + PropertyManager.get_instance().get_property( DataConstants.HEADER_FILE_OUTPUT_PROPERTY, True)
        
        self.source = PropertyManager.get_instance().get_property( DataConstants.HEADER_SOURCE_PROPERTY, True)
        self.type_id = PropertyManager.get_instance().get_property( DataConstants.HEADER_TYPE_ID_PROPERTY,True)
        
        # Method calling
        HeaderParser.change_header(self.path_file_input, self.path_file_output, source=int(self.source), type_id=int(self.type_id))
        
    
    # split_dataset
    # -----------------
    #         
    # Division of sequences file in many files each containing the sequence of one protein
    # This step is necessary if you want analyze the proteins sequences with anchor and Iupred
    #
    # Output: As many fasta file as the proteins are in the fasta file input will be created in the path_output
    #     
    def split_dataset(self):
            
            
        Logger.get_instance().info( " Division of Dataset in many fasta file each containing one protein sequence")
        
        self.path_home = Constants.PATH_HOME
            
        self.split_path_input = self.path_home + PropertyManager.get_instance().get_property( DataConstants.SPLIT_PATH_INPUT_PROPERTY, True)
        self.split_file_fasta = self.split_path_input + PropertyManager.get_instance().get_property( DataConstants.SPLIT_DATASET_PROPERTY, True)
        self.split_path_output = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.SPLIT_PATH_OUTPUT_PROPERTY, True)
        self.split_start_index = PropertyManager.get_instance().get_property( DataConstants.SPLIT_START_HEADER_PROPERTY, True)
        self.split_end_index = PropertyManager.get_instance().get_property( DataConstants.SPLIT_END_HEADER_PROPERTY, True)
        
        SplitSeq.split_seq( self.split_file_fasta, self.split_path_output, int(self.split_start_index), int(self.split_end_index))
        
            
        Logger.get_instance().info( " The dataset has been split in many fasta files ")
            
            
    # anchor_analysis
    # -----------------
    # 
    # This method performs the Anchor Analysis 
    # The global anchor analysis method is called
    #
    # Output: The Anchor output files are stored in anchor_path_output (PBR_proteinName.txt,FMotifs_proteinName.txt, Pred_proteinName.txt)
    #
    def anchor_analysis(self):
        
        
        
        Timer.get_instance().step(" Start of Anchor analysis...")
        
        self.tool_path_input = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.TOOL_PATH_INPUT_PROPERTY, True)
        self.anchor_path_output = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.ANCHOR_PATH_OUTPUT_PROPERTY,True)
        self.motif_list_path = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.ANCHOR_MOTIF_PATH_PROPERTY, True)
        
        Anchor.global_anchor_analysis(self.motif_list_path,
                                      self.tool_path_input,
                                      self.anchor_path_output)
        
        Timer.get_instance().step(" End of Anchor analysis")
        
        
        
    # iupred_analysis
    # ----------------
    # 
    # This method performs Iupred anaysis 
    # The global iupred analysis is called
    #
    # Output: The Iupred output is stored in iupred_path_ouptut (IUPred_ProteName.txt)
    #  
    def iupred_analysis(self):
        
        Timer.get_instance().step(" Start of Iupred analysis...")
        
        
        self.tool_path_input = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.TOOL_PATH_INPUT_PROPERTY, True)
        self.iupred_path_output = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.IUPRED_PATH_OUTPUT_PROPERTY,True)
        
            
        Iupred.global_iupred_analysis(self.tool_path_input,
                                      self.iupred_path_output)
        
        Timer.get_instance().step(" End of Iupred analysis")
        
    
    
    # analysis_tools_output
    # ----------------------------
    #
    # This method carries out the Analysis of Iupred and Anchor output
    # producing some files containing the information that will can be drawn 
    # 
    #
    def analysis_tools_output(self):
        
        Timer.get_instance().step(" Start of tools analysis.. ")
        
        
        self.path_home = Constants.PATH_HOME
        self.input_path_iupred = self.path_home + PropertyManager.get_instance().get_property( DataConstants.ANALYSIS_INPUT_PATH_IUPRED_PROPERTY, True)
        self.output_path_analysis = self.path_home + PropertyManager.get_instance().get_property( DataConstants.ANALYSIS_OUTPUT_PATH_TOOLS_PROPERTY, True)
        
        self.threshold_1 =  PropertyManager.get_instance().get_property( DataConstants.ANALYSIS_THRESHOLD_1_PROPERTY, True)
        self.threshold_2 =  PropertyManager.get_instance().get_property( DataConstants.ANALYSIS_THRESHOLD_2_PROPERTY, True)
        self.number_aa_iupred =  PropertyManager.get_instance().get_property( DataConstants.ANALYSIS_AMINOACID_NUMBER_IUPRED_PROPERTY, True)
        self.dataset_type = PropertyManager.get_instance().get_property( DataConstants.ANALYSIS_DATASET_TYPE_PROPERTY, True)
        
        Iupred.make_iupred_file(self.input_path_iupred, self.output_path_analysis, float(self.threshold_1), float(self.threshold_2), int(self.number_aa_iupred), self.dataset_type)
        
        self.input_path_anchor = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.ANALYSIS_INPUT_PATH_ANCHOR_PROPERTY, True)
        self.num_aa_anchor = PropertyManager.get_instance().get_property( DataConstants.ANALYSIS_AMINOACID_NUMBER_ANCHOR_PROPERTY, True)
             
        Anchor.make_anchor_file(self.input_path_anchor, self.output_path_analysis, int(self.num_aa_anchor),self.dataset_type )
                
        
        Timer.get_instance().step(" End of tools analysis")
        
    
    # disordpbind_analysis
    # --------------------------
    # 
    # This method carries out the Analysis of DisoRDPbind output 
    # DisoRDPbind tool can not installed in local then the output file must be provided to run this analysis 
    #
    # 
    
    def disordpbind_analysis(self):
        
        Timer.get_instance().step(" Start of DisoRDPbind output analysis.. ")
        
        self.path_home = Constants.PATH_HOME
        self.input_file = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DISO_INPUT_FILE_PROPERTY, True)
        self.ouput_path = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DISO_OUTPUT_FOLDER_PROPERTY, True)
        self.binding_partner = PropertyManager.get_instance().get_property( DataConstants.DISO_BINDING_PARTNER_PROPERTY, True)
        self.num_aa_diso = PropertyManager.get_instance().get_property( DataConstants.DISO_NUM_AA_PROPERTY, True)
        self.dataset_type = PropertyManager.get_instance().get_property( DataConstants.DISO_DATASET_TYPE_PROPERTY, True)
        
        
        
        DisoRDPbind.make_disordp_file(self.input_file, self.ouput_path, int(self.binding_partner), int(self.num_aa_diso), self.dataset_type)
        
        Timer.get_instance().step(" End of DisoRDPbind output analysis")
        
    
    # particular_analysis
    # -------------------------------
    # 
    # This method has been created in specific way for RBP dataset because the RBP dataset has been classified in different way:
    # in according to type of RNA and to class of domain
    # Nevertheless it can be applied also to other dataset providing lists with specific protein classification  
    # 
    # This method takes as input the Analysis file of Iupred, Anchor and DisoRDPbind and classify the result 
    # in according to proteins lists provided as input.
    #
    # 
    
    def particular_analysis(self):
        
        
        Timer.get_instance().step(" Start of tools analysis for specific protein ")
        
        self.path_home = Constants.PATH_HOME
        self.path_input_anchor_file = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_INPUT_ANCHOR_FILE_PROPERTY, True)
        self.path_input_iupred_file = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_INPUT_IUPRED_FILE_PROPERTY, True)
        self.path_input_disordp_file = self.path_home + PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_INPUT_DISORDP_FILE_PROPERTY, True)
        self.path_input_reg_anchor = self.path_home + PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_INPUT_REG_ANCHOR_FILE_PROPERTY, True)
        self.path_input_reg_iupred_1 = self.path_home + PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_INPUT_REG_IUPRED_1_FILE_PROPERTY, True)
        self.path_input_reg_iupred_2 = self.path_home + PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_INPUT_REG_IUPRED_2_FILE_PROPERTY, True)
        self.path_input_reg_diso = self.path_home + PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_INPUT_REG_DISO_FILE_PROPERTY, True)
        self.input_files = self.path_home + PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_INPUT_DIR_FILE_PROPERTY, True)
        self.list_namefiles = PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_LIST_NAMEFILE_PROPERTY, True)
        self.path_output_dir = self.path_home + PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_OUTPUT_DIR_PROPERTY, True)
        self.path_output_dir_diso = self.path_home + PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_OUTPUT_DIR_DISO_PROPERTY, True)
        
        # This parameter represents the column of protein id in the classification files
        #
        # In Domain Class  files the column of protein id is the 2 ( that is 1 for python)
        # In RNA target files the column of protein id is the 1 (that is 0 for python)
        #
        
        #self.protein_column_rna =  PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_PROTEIN_LIST_COLUMN_RNA_PROPERTY, True)
        self.protein_column_class =  PropertyManager.get_instance().get_property( DataConstants.SPECIFIC_PROTEIN_LIST_COLUMN_CLASS_PROPERTY, True)
        
        # region file
        anchor_table = FileParser.make_table(self.path_input_reg_anchor, skip=1)
        iupred_table_1 = FileParser.make_table(self.path_input_reg_iupred_1, skip=1)
        iupred_table_2 = FileParser.make_table(self.path_input_reg_iupred_2, skip=1)
        disordp_table = FileParser.make_table(self.path_input_reg_diso, skip=1)
        
        # table file (fraction)
        anchor_t = FileParser.make_table(self.path_input_anchor_file, skip=1)
        iupred_t = FileParser.make_table(self.path_input_iupred_file, skip=1)
        disordp_t = FileParser.make_table(self.path_input_disordp_file, skip=1)
        
        list_filenames = self.list_namefiles.split(',')
        
        for filename in list_filenames:
            feature = filename.split('.')[0]
            table_domain = FileParser.make_table(self.input_files + str(filename))
            list_prot = TableWrapper.get_column(table_domain,int(self.protein_column_class))
            prot_id_anchor = TableWrapper.get_column(anchor_table, 0)
            prot_id_iupred_1 = TableWrapper.get_column(iupred_table_1, 0)
            prot_id_iupred_2 = TableWrapper.get_column(iupred_table_2, 0)
            prot_id_disordp = TableWrapper.get_column(disordp_table, 0)
            
            prot_id_anchor_t = TableWrapper.get_column(anchor_t, 0)
            prot_id_iupred_t = TableWrapper.get_column(iupred_t, 0)
            prot_id_disordp_t = TableWrapper.get_column(disordp_t, 0)
            
            # region file
            new_table_anchor = [line for n, line in enumerate(anchor_table) if prot_id_anchor[n] in list_prot]
            new_table_iupred_1 = [line for n, line in enumerate(iupred_table_1) if prot_id_iupred_1[n] in list_prot]
            new_table_iupred_2 = [line for n, line in enumerate(iupred_table_2) if prot_id_iupred_2[n] in list_prot]
            new_table_disordp = [line for n, line in enumerate(disordp_table) if prot_id_disordp[n] in list_prot]
            anchor_output_file_path = self.path_output_dir + feature + '_AnchorRegion.txt'
            iupred1_output_file_path = self.path_output_dir + feature + '_IUPredRegion_0.4.txt'
            iupred2_output_file_path = self.path_output_dir + feature + '_IUPredRegion_0.5.txt'
            disordp_output_file_path = self.path_output_dir_diso + feature + '_DisoRDPRegion.txt'
            
            # Table file (fraction)
            new_table_a = [line for n, line in enumerate(anchor_t) if prot_id_anchor_t[n] in list_prot]
            new_table_i = [line for n, line in enumerate(iupred_t) if prot_id_iupred_t[n] in list_prot]
            new_table_d = [line for n, line in enumerate(disordp_t) if prot_id_disordp_t[n] in list_prot]
            anchor_output_table = self.path_output_dir + feature + '_AnchorTable.txt'
            iupred_output_table = self.path_output_dir + feature + '_IUPredTable_0.4_0.5.txt'
            disordp_output_table = self.path_output_dir_diso + feature + '_DisoRDPTable.txt'
            
            
            # file writing
            
            # Region file
            FileWriter.write_table(anchor_output_file_path, new_table_anchor)
            FileWriter.write_table(iupred1_output_file_path, new_table_iupred_1)
            FileWriter.write_table(iupred2_output_file_path, new_table_iupred_2)
            FileWriter.write_table(disordp_output_file_path, new_table_disordp)
            
            # Table file
            FileWriter.write_table(anchor_output_table, new_table_a)
            FileWriter.write_table(iupred_output_table, new_table_i)
            FileWriter.write_table(disordp_output_table, new_table_d)
        
        Timer.get_instance().step(" End of tools analysis for specific protein ")

    
    # whole_procedure
    # ----------------
    # 
    # This method allows to select just some of previously methods
    #

    @staticmethod
    def whole_procedure():
            
        
        # start chrono
        Timer.get_instance().start_chrono()
        Logger.get_instance().info("Start of Disorder Analysis.....\n ")
           
        disorder = DisorderAnalysis()
        #disorder.change_header()
        #disorder.split_dataset()
        #disorder.anchor_analysis()
        #disorder.iupred_analysis()
        #disorder.analysis_tools_output()
        #disorder.disordpbind_analysis()
        #disorder.particular_analysis()
            
        Timer.get_instance().stop_chrono(' End of Disorder Analysis')
        
     


if __name__ == '__main__':
    
    OptionManager.get_instance().initialize()

    PropertyManager.get_instance().read_properties( OptionManager.get_instance().get_option( OptionConstants.OPTION_DISORDER_ANALYSIS_PATH, True))
    
    DisorderAnalysis.whole_procedure()
    
    