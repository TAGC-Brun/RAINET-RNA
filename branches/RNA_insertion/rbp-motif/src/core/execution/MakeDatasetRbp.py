# ===========================================================================
# This class contains the steps for creation of RBP protein sequences dataset
# ===========================================================================


# Starting file
# Dataset 1 
# Dataset 2

from core.util import Constants
from core.util.property.PropertyManager import PropertyManager
from core.util.option import OptionConstants
from core.util.option.OptionManager import OptionManager
from core.util.log.Logger import Logger
from core.data import DataConstants
from core.util.dataset.InfoDataset import InfoDataset
from core.util.time.Timer import Timer
from core.util.databases.Ensembl import Ensembl
from core.util.file.FileUtils import FileUtils
from core.util.parsing.TableWrapper import TableWrapper
from core.util.parsing.FileParser import FileParser
from core.util.file.FileWriter import FileWriter
from core.util.format.InfoFasta import InfoFasta
import pickle
from core.util.format.LengthSeq import LengthSeq
import random
from core.util.format.SeqToFasta import SeqToFasta
from core.util.file.RemoveFile import RemoveFile


class MakeDatasetRbp():
        
    # comparison_dataset
    # --------------------
    #
    # This method compares two datasets producing two list file containing the common and different gene 
    # (The different item are the item in dataset b and not in dataset a)
    #
    # 
    def comparison_dataset(self):
        
        
        # Create Logger instance to see the start of comparison between two datasets 
        Logger.get_instance().info( " Start of comparison datasets:..." )
        
        
        # Definition of InfoDataset arguments
        
        self.dataset_input_path = PropertyManager.get_instance().get_property( DataConstants.DATASET_INPUT_PATH_PROPERTY, True)
        self.dataset_1_file = PropertyManager.get_instance().get_property( DataConstants.DATASET_1_FILE_PROPERTY, True)
        self.dataset_2_file =  PropertyManager.get_instance().get_property( DataConstants.DATASET_2_FILE_PROPERTY, True)
        self.dataset_1_index_col = PropertyManager.get_instance().get_property( DataConstants.DATASET_1_INDEX_COL_PROPERTY, True)
        self.dataset_2_index_col = PropertyManager.get_instance().get_property( DataConstants.DATASET_2_INDEX_COL_PROPERTY, True)
        self.dataset_output = PropertyManager.get_instance().get_property( DataConstants.DATASET_OUTPUT_PROPERTY, True)
        self.dataset_1_length = PropertyManager.get_instance().get_property( DataConstants.DATASET_1_LENGTH_PROPERTY, True)
        self.dataset_2_length = PropertyManager.get_instance().get_property( DataConstants.DATASET_2_LENGTH_PROPERTY, True)
        
        self.index_col = (int(self.dataset_1_index_col) , int(self.dataset_2_index_col))
        self.dataset_length = (int(self.dataset_1_length), int(self.dataset_2_length))
        
        self.path_home = Constants.PATH_HOME
        InfoDataset.global_analysis_dataset(self.path_home + self.dataset_input_path, 
                                           (self.dataset_1_file,
                                           self.dataset_2_file) ,
                                           self.index_col,
                                           self.path_home + self.dataset_output,
                                           length=self.dataset_length)
        
        Logger.get_instance().info( " The comparison of datasets is completed : \
two file  with the common and difference \
items has been generated in \n\n " + self.dataset_output )
        
        
    # creation_list
    # ----------------
    #
    # Creation of the gene and protein list of dataset 1 in order to pass it to Ensembl method 
    #    
    def creation_list(self):
        
        
        
        Logger.get_instance().info( " Creation of gene and protein list \n " )
        
        
        # Creation of two file containing respectively the genes and protein of dataset 1
        
        self.dataset_input_path = PropertyManager.get_instance().get_property( DataConstants.DATASET_INPUT_PATH_PROPERTY, True)
        self.file_dataset_1 = PropertyManager.get_instance().get_property( DataConstants.DATASET_1_FILE_PROPERTY, True)
        self.gene_index_col = PropertyManager.get_instance().get_property( DataConstants.LIST_GENE_INDEX_COL_PROPERTY, True)
        self.protein_index_col = PropertyManager.get_instance().get_property( DataConstants.LIST_PROTEIN_INDEX_COL_PROPERTY, True)
        
        
        self.list_gene_dataset_1 = PropertyManager.get_instance().get_property( DataConstants.LIST_FILE_GENE_DATASET_1_PROPERTY, True)
        self.list_protein_dataset_1 = PropertyManager.get_instance().get_property( DataConstants.LIST_FILE_PROTEIN_DATASET_1_PROPERTY, True)
        
        self.path_gene_dataset_1 = Constants.PATH_HOME + self.dataset_input_path + self.list_gene_dataset_1
        self.path_protein_dataset_1 = Constants.PATH_HOME + self.dataset_input_path + self.list_protein_dataset_1
        
        self.path_home = Constants.PATH_HOME
        self.path_dataset_1 = self.path_home + self.dataset_input_path + self.file_dataset_1
        
        dataset_1 = FileParser.make_table(self.path_dataset_1)
        
        gene_dataset_1 = TableWrapper.get_column(dataset_1, int(self.gene_index_col), start=1)
        protein_dataset_1 = TableWrapper.get_column(dataset_1, int(self.protein_index_col), start=1)
        
        FileWriter.write_table(self.path_gene_dataset_1, gene_dataset_1)
        FileWriter.write_table(self.path_protein_dataset_1, protein_dataset_1)
        
        
        Logger.get_instance().info( " The genes and proteins file of dataset 1 have been created \
in the following path  \n\n " + self.dataset_input_path )
        
        
    # connection_to_ensembl
    # ----------------------
    #  
    # Section reserved to extraction of sequences from Ensembl database    
    # DATASET 1:
    # This dataset contains the gene and protein ids then you can download directly the protein sequence
    # 
    # DATASET 2:
    # For this dataset are download only the sequences for the gene that not are also in Dataset 1 
    # That is the Differ items obtained from comparison_dataset method
    #
    def connection_to_ensembl(self):
          
        Logger.get_instance().info( " Connection to Ensembl: Starting...\n" )
        
        # DATASET 1
        # =============================================
        
        # Collection of sequences for dataset 1 
        Logger.get_instance().info( " Dataset 1 : Extraction of sequences...\n" )
        
        # Timer step
        Timer.get_instance().step(" Start of sequences extraction \n")
    
        
        # Definition of Ensembl list_get_seq arguments 
        
        self.path_home = Constants.PATH_HOME
        self.list_path = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.DATASET_INPUT_PATH_PROPERTY, True)
        
        self.gene_list_1 = PropertyManager.get_instance().get_property( DataConstants.LIST_FILE_GENE_DATASET_1_PROPERTY, True)
        self.protein_list = PropertyManager.get_instance().get_property( DataConstants.LIST_FILE_PROTEIN_DATASET_1_PROPERTY, True)
        
        self.ensembl_gene_list_1_path = self.list_path + self.gene_list_1
        self.ensembl_protein_list_1_path = self.list_path + self.protein_list
        self.type_query1_ensembl = PropertyManager.get_instance().get_property( DataConstants.ENSEMBL_TYPE_QUERY_DATASET_1_PROPERTY, True)
        
        
        self.ensembl_path_output = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.ENSEMBL_OUTPUT_PATH_SEQUENCE_PROPERTY, True)
        self.ensembl_output_dataset1 =  self.ensembl_path_output +  PropertyManager.get_instance().get_property( DataConstants.ENSEMBL_FILE_SEQUENCES_1_PROPERTY)



        # Calling Ensembl.list_get_seq
        Ensembl.list_get_seq(self.ensembl_gene_list_1_path, int(self.type_query1_ensembl), path_protein_list=self.ensembl_protein_list_1_path, path_output=self.ensembl_output_dataset1)
        
        
        
        # Timer step
        Timer.get_instance().step(" End of Dataset 1 Sequences Extraction\n")
        

        Logger.get_instance().info( " Extraction of sequences for the dataset 1 has been completed \n\n" )
        
        
        # END DATASET 1
        # =====================================================
        
        
        
        
        
        # DATASET 2
        # =====================================================
        

        # Collection of sequences for dataset 2 
        Logger.get_instance().info( " Dataset 2 : Extraction of sequences ....\n" )
        
    
        
        # Definition of Ensembl list_get_seq arguments 
        
        
        self.ensembl_input_list_2 = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.DATASET_OUTPUT_PROPERTY, True)
        self.gene_list_2 = Constants.FILE_DIFF
        self.ensembl_gene_list_2_path = self.ensembl_input_list_2 + self.gene_list_2
        self.type_query2_ensembl = PropertyManager.get_instance().get_property( DataConstants.ENSEMBL_TYPE_QUERY_DATASET_2_PROPERTY, True)
        
        self.ensembl_path_output = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.ENSEMBL_OUTPUT_PATH_SEQUENCE_PROPERTY, True)
        self.ensembl_output_dataset2 = self.ensembl_path_output + PropertyManager.get_instance().get_property( DataConstants.ENSEMBL_FILE_SEQUENCES_2_PROPERTY, True)
       
        
        # Calling Ensembl.list_get_seq
        Ensembl.list_get_seq(self.ensembl_gene_list_2_path, int(self.type_query2_ensembl), path_protein_list=None, path_output = self.ensembl_output_dataset2)
        
        
        # Timer step
        Timer.get_instance().step(" End of Dataset 2 Sequences Extraction\n")
        

        Logger.get_instance().info( " Extraction of sequences for the dataset 2 has been completed \n\n" )
        
        # END DATASET 2
        # =====================================================
        
        
        
        Logger.get_instance().info( " The sequences file of dataset 1 and the novel gene in dataset 2 \
have been created in the following path  \n" + self.ensembl_path_output)
          
    # dictionary_identifier
    # ----------------------
    # 
    # Creation of a dictionary: gene = [ isoform1, isoform2,...isoformN] 
    # This method is necessary to select the longest isoform for each gene    
    # 
    def dictionary_identifier(self):
        
        
        Logger.get_instance().info( " Creation of a dictionary for novel gene of dataset 2\
The dictionary structure is : \n \
{gene = [ isoform1, isoform2,...isoformN]}" )
        
        self.ensembl_path_output = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.ENSEMBL_OUTPUT_PATH_SEQUENCE_PROPERTY, True)
        self.ensembl_output_dataset2 = self.ensembl_path_output + PropertyManager.get_instance().get_property( DataConstants.ENSEMBL_FILE_SEQUENCES_2_PROPERTY, True)
       

        self.dictionary_output = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.DICTIONARY_PATH_OUTPUT_PROPERTY, True)
        self.dictionary_namefile = self.dictionary_output + PropertyManager.get_instance().get_property( DataConstants.DICTIONARY_NAME_FILE_PROPERTY, True) 
        
        dict_identifier = InfoFasta.make_dictionary(self.ensembl_output_dataset2)
        
        file_dict = FileUtils.open_text_w(self.dictionary_namefile)
        
        pickle.dump(dict_identifier, file_dict)
        
        
        Logger.get_instance().info( " The creation of a dictionary for novel gene in dataset 2 is completed \n\n")
        
    # longest_sequence
    # ----------------
    # 
    # this method find the longest isoform for each gene (see the principals of LengthSeq method)
    # 
    def longest_sequence(self):
        
        
        Logger.get_instance().info( " Start of the selection of longest sequences of novel dataset \n")
        
        
        
        # Definition of arguments
        self.path_sequences = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.LONGEST_PATH_SEQUENCE_PROPERTY, True)
        self.file_sequences =  self.path_sequences + PropertyManager.get_instance().get_property( DataConstants.LONGEST_PROT_FILE_SEQUENCES_2_PROPERTY, True)
        self.path_dictionary_identifier = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.LONGEST_PATH_DICTIONARY_PROPERTY, True)
        self.file_dictionary = self.path_dictionary_identifier + PropertyManager.get_instance().get_property( DataConstants.LONGEST_DICTIONARY_NAME_FILE_PROPERTY, True)
        
        self.path_output_longest = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.LONGEST_PATH_OUTPUT_PROPERTY, True)
        
        self.path_file_longest = self.path_output_longest + PropertyManager.get_instance().get_property( DataConstants.LONGEST_FILE_PROPERTY, True)
        self.path_file_isoform = self.path_output_longest + PropertyManager.get_instance().get_property( DataConstants.ISOFORM_FILE_PROPERTY, True)
        
        
        # Extraction the longest sequences from dataset 2 sequences (isoforms)  
        LengthSeq.longest_seq(self.file_sequences, self.file_dictionary, self.path_file_longest, self.path_file_isoform)
        
        
        # Timer step
        Timer.get_instance().step(" End of selection of the longest sequences in dataset 2  \n")
        
        
        Logger.get_instance().info( " End of selection of the longest sequences: \n \
two file have been generated one with longest sequences and the other one containing the isoform with same length  ")
    
    
    
         
    # isoform_sequences 
    # --------------------------
    #
    # like explained in LengthSeq method a file containing isoform with same length is created
    #
    # In this section occurs the Random Selection of isoforms sequences
    # 
    def isoform_sequences(self):
        
        
        Logger.get_instance().info( " Starting the random selection of isoforms with same length \n")
        Logger.get_instance().info( " The following headers are the proteins randomly selected \n")
        
        self.path_output_longest = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.LONGEST_PATH_OUTPUT_PROPERTY, True)
        
        
        self.path_file_isoform = self.path_output_longest + PropertyManager.get_instance().get_property( DataConstants.ISOFORM_FILE_PROPERTY, True)
        self.path_file_selected_isoform = self.path_output_longest + PropertyManager.get_instance().get_property( DataConstants.RANDOM_ISOFORM_SEQ_PROPERTY, True)
        
        # The headers of a Isoform fasta file are taken by InfoFasta class
        # You make sure that the arg text is equal to False because the input object is a file and not a list
                                                                                                              
                                                                                                
        self.headers = InfoFasta.get_header(self.path_file_isoform,  text=False)
        
        # Extraction of genes form headers line
        # This vector contains double gene because the file contains some isoform of the same gene
        
        gene_isoform = []
        for header in self.headers:
            gene = header[1:16]
            gene_isoform.append(gene)
        
        # gene set creation    
        unique_gene = set(gene_isoform)
        
        # This for loop flows on the unique gene 
        #
        random_header = []
        old_num_isoform = 0
        for gene in unique_gene:
            # For each gene counts how many isoform has 
            num_isoform = gene_isoform.count(gene)
            item = range(0,num_isoform)
            # Select one isoform randomly
            sel = random.choice(item)
            # The header selected randomly are stored in array
            random_header.append(self.headers[old_num_isoform : old_num_isoform + num_isoform][sel])
            old_num_isoform = old_num_isoform + num_isoform
            
            
        self.file_random_seq = FileUtils.open_text_a(self.path_file_selected_isoform)
        
        # The sequences corresponding to header selected are extracted from isoform file   
            
        for header in random_header:
            Logger.get_instance().info('Header selected : ' + header)
            identifier = header[33:48]
            sequence = InfoFasta.get_seq(self.path_file_isoform, identifier)
            fasta_seq = SeqToFasta.give_fasta(header, sequence)
            self.file_random_seq.write(fasta_seq)
            
        
        Logger.get_instance().info( " End of selection random sequences \n ")
        
       
        
    # merger_sequences
    # --------------------
    #
    # This method merges:
    #    - the Longest sequences and the random isoform (dataset 2)
    #    - FinalSeqDataset2 + SequencesDataset1
    # 
    # The output is the final fasta file of dataset that can be gived as input to DisorderAnalysis
    #    


    def merger_sequences(self):
            
        Logger.get_instance().info( " Union of the longest sequences and the random selected isoform ")
            
        # Input variables to merge the longest Novel sequences and random selected isoform of dataset 2
        
        self.path_input_longest = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.FUSION_PATH_INPUT_PROPERTY, True)
        self.path_file_longest = self.path_input_longest + PropertyManager.get_instance().get_property( DataConstants.LONGEST_FILE_PROPERTY, True)
        self.path_file_isoform = self.path_input_longest + PropertyManager.get_instance().get_property( DataConstants.SELECTED_ISOFORM_FILE_PROPERTY, True)
        
        self.path_output_seq = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.FUSION_PATH_OUTPUT_PROPERTY, True)
        self.path_file_seq_dataset_2 = self.path_output_seq +  PropertyManager.get_instance().get_property( DataConstants.FUSION_FILE_SEQ_DATASET_2_PROPERTY, True)
        
        
        FileParser.merge_file(self.path_file_longest, self.path_file_isoform, self.path_file_seq_dataset_2)
        
                
        # Input variables to merge the sequences datasets (Novel_JProteomics and NatRevGenetics)
        
        self.path_input_seq_dataset1 = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.FUSION_PATH_INPUT_DATASET_1_PROPERTY, True)
        self.path_file_dataset1 = self.path_input_seq_dataset1 +  PropertyManager.get_instance().get_property( DataConstants.FUSION_FILE_DATASET_1_PROPERTY, True)
        
        self.path_file_dataset12 = self.path_output_seq + PropertyManager.get_instance().get_property( DataConstants.FUSION_DATASET_12_PROPERTY, True)
            
            
            
        Logger.get_instance().info( " Union of sequences respectively of dataset 1 and the novel dataset 2 proteins \n ")
            
        FileParser.merge_file(self.path_file_dataset1, self.path_file_seq_dataset_2, self.path_file_dataset12)
            
        
        Logger.get_instance().info( " The New RBP Dataset has been created\n ")
            
            
        # This part checks if there are pseudo genes inside dataset 2
        # Make the comparison between the original genes gived as an input and the genes obtained after 
        # connection to Ensembl 
        # This check allows to find genes that are not anymore available or that are pseudogenes
         
        Logger.get_instance().info( " Comparison between original genes and Ensembl output ")
        self.path_home = Constants.PATH_HOME
        self.path_input_original_file = self.dataset_output = PropertyManager.get_instance().get_property( DataConstants.DATASET_OUTPUT_PROPERTY, True)
        self.original_file = self.path_home + self.path_input_original_file + Constants.FILE_DIFF
        
        original_genes = FileParser.read_file(self.original_file)
                
        
        self.path_output_seq = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.FUSION_PATH_OUTPUT_PROPERTY, True)
        self.path_file_seq_dataset_2 = self.path_output_seq +  PropertyManager.get_instance().get_property( DataConstants.FUSION_FILE_SEQ_DATASET_2_PROPERTY, True)
        
        final_headers = InfoFasta.get_header(self.path_file_seq_dataset_2)
        final_genes = [item[1:16] for item in final_headers]
        
        out_comparison = InfoDataset.comparison_dataset(original_genes, final_genes, header=False)
        
        genes = '\n'.join(out_comparison[1])
        
        Logger.get_instance().info(" The genes lost during the request to Ensembl are : \n" + genes)
        

       
        
        
    
    
    # method that delete the append file in order to run again the whole workflow
    def delet_append_file(self):
        
        
        self.del_ensembl_input = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.DEL_ENSEMBL_PATH_PROPERTY, True)
        self.del_ensembl_file1 =  self.del_ensembl_input + PropertyManager.get_instance().get_property( DataConstants.DEL_ENSEMBL_FILE1_PROPERTY, True)
        self.del_ensembl_file2 =  self.del_ensembl_input + PropertyManager.get_instance().get_property( DataConstants.DEL_ENSEMBL_FILE2_PROPERTY, True)
        
        self.del_longest_input = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.DEL_LONGEST_PATH_PROPERTY, True)
        self.del_longest_file =  self.del_longest_input + PropertyManager.get_instance().get_property( DataConstants.DEL_LONGEST_FILE_PROPERTY, True)
        self.del_isoform_file =  self.del_longest_input + PropertyManager.get_instance().get_property( DataConstants.DEL_ISOFORM_FILE_PROPERTY, True)
        self.del_random_isoform_file = self.del_longest_input + PropertyManager.get_instance().get_property( DataConstants.DEL_RANDOM_ISOFORM_FILE_PROPERTY, True)
        
        self.del_fusion_path = Constants.PATH_HOME + PropertyManager.get_instance().get_property( DataConstants.DEL_FUSION_PATH_PROPERTY, True)
        self.del_fusion_file_longest = self.del_fusion_path + PropertyManager.get_instance().get_property( DataConstants.DEL_FUSION_DATASET_LONGEST_PROPERTY, True)
        self.del_fusion_file_dataset12 = self.del_fusion_path + PropertyManager.get_instance().get_property( DataConstants.DEL_FUSION_DATASET12_PROPERTY, True)
        
        files = [self.del_ensembl_file1, 
                 self.del_ensembl_file2,
                 self.del_longest_file,
                 self.del_isoform_file,
                 self.del_random_isoform_file,
                 self.del_fusion_file_longest,
                 self.del_fusion_file_dataset12]
        
        for namefile in files:
            RemoveFile.delete_file(namefile)
        

        
        
    # whole_procedure
    # ----------------
    # 
    # This method allows to select just some of previously methods
    #
        
    @staticmethod
    def whole_procedure():
            
        
        # start chrono
        Timer.get_instance().start_chrono()
        Logger.get_instance().info("Start of the creation of RBP dataset.....\n ")
           
        M = MakeDatasetRbp()
        
        #M.delet_append_file()
        
        #M.comparison_dataset()
        #M.creation_list()
        #M.connection_to_ensembl()
        #M.dictionary_identifier()
        #M.longest_sequence()
        #M.isoform_sequences()
        #M.merger_sequences()
        M.split_dataset()
            
        Timer.get_instance().stop_chrono(' End of the creation of RBP dataset')
        
        
        
if __name__ == '__main__':
    
    OptionManager.get_instance().initialize()
    
    # Retrieve the MakeDatasetRbp properties
    PropertyManager.get_instance().read_properties( OptionManager.get_instance().get_option( OptionConstants.OPTION_MAKEDATASETRBP_PROPERTIES_PATH, True))

    
    MakeDatasetRbp.whole_procedure()

            
             
             
             
             
                               
            
        
        
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        