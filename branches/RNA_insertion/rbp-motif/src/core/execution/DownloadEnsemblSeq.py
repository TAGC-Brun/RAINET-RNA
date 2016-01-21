# ================================================================================
# This class takes the longest sequences from Ensembl giving a genes list as input
# ================================================================================




from core.util import Constants
from core.util.property.PropertyManager import PropertyManager
from core.data import DataConstants
from core.util.option import OptionConstants
from core.util.option.OptionManager import OptionManager
from core.util.time.Timer import Timer
from core.util.databases.Ensembl import Ensembl
from core.util.format.InfoFasta import InfoFasta
from core.util.file.FileUtils import FileUtils
from core.util.format.LengthSeq import LengthSeq
import random
from core.util.format.SeqToFasta import SeqToFasta
from core.util.parsing.FileParser import FileParser

from core.util.log.Logger import Logger
import pickle




class DownloadEnsemblSeq():
    
    
    # download_product_gene_seq
    # ----------------------------
    #
    # Method reserved to extraction of sequences from Ensembl database  
    # It takes as input a list of gene
    # 
    def download_product_gene_seq(self):
        
        
        Logger.get_instance().info( " Start of sequences download from Ensembl..")
        
        self.path_home = Constants.PATH_HOME
        self.gene_list_input = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_ENSEMBL_FILE_INPUT_LIST_PROPERTY, True)
        self.ensembl_seq_output = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_ENSEMBL_FILE_OUPUT_SEQ_PROPERTY, True)
        self.type_query = PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_ENSEMBL_TYPE_QUERY_PROPERTY, True)
        
        
        Ensembl.list_get_seq(self.gene_list_input, int(self.type_query), path_output=self.ensembl_seq_output)


        Logger.get_instance().info( " End of sequences download from Ensembl..")
        
    
    # make_dictionary
    # ----------------------
    #
    # Creation of a dictionary: gene = [ isoform1, isoform2,...isoformN] 
    # This method has to use only for the sequences downloaded from ensembl if you have a list of sequences 
    # you can skip this workflow and go directly to Disorder Analysis workflow
    #
    # Note: If the sequences are download manually and not automatically from BioMart you make sure that the sequences 
    # haven't the asterisk
    #
    def make_dictionary(self):
        
        Logger.get_instance().info( " Creation of a dictionary for novel gene of dataset 2\
The dictionary structure is : \n \
{gene = [ isoform1, isoform2,...isoformN]}" )
        
        
        self.path_home =  Constants.PATH_HOME
        self.path_input_file = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_DICTIONARY_INPUT_FILE_PROPERTY, True)
        
        self.dictionary_output_path = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_DICTIONARY_OUTPUT_PATH_PROPERTY, True)
        self.output_file_path = self.dictionary_output_path + PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_DICTIONARY_FILE_OUTPUT_PROPERTY, True)
        
        
        dict_identifier = InfoFasta.make_dictionary(self.path_input_file)
        
        self.dict_file = FileUtils.open_text_w(self.output_file_path)
        
        pickle.dump(dict_identifier, self.dict_file)       
    
        Logger.get_instance().info( " The creation of a dictionary is completed \n\n")
        
        
    # get_longest_seq
    # -----------------
    #
    # this method find the longest isoform for each gene (see the principals of LengthSeq method)
    # 
    #
    def get_longest_seq(self):
        
        
        Logger.get_instance().info( " Start of the selection of longest sequences of dataset \n")
        
        
        self.path_home = Constants.PATH_HOME
        
        self.file_sequences = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_LONGEST_SEQ_INPUT_FILE_PROPERTY, True)
        self.output_path = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_LONGEST_SEQ_OUTPUT_PATH_PROPERTY, True)
        
        self.dictionary_file = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_LONGEST_DICTIONARY_PROPERTY, True)
        
        self.longest_file = self.output_path +  PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_LONGEST_SEQ_FILE_PROPERTY, True)
        self.isoform_file = self.output_path +  PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_LONGEST_ISOFORM_FILE, True)
        
        
        # Extraction the longest sequences from dataset sequences (isoforms)  
        self.file_seq = open(self.file_sequences)
        self.seq_obj = self.file_seq.readlines()
        
        LengthSeq.longest_seq(self.seq_obj, self.dictionary_file, self.longest_file, self.isoform_file, type_obj='list')
        
        
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
        
        
        self.path_home = Constants.PATH_HOME
        self.path_file_isoform = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_ISOFORM_FILE_PATH_PROPERTY, True)
        
        self.path_random_file = self.path_home +   PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_RANDOM_FILE_PATH_PROPERTY, True)
        
        # The headers of a Isoform fasta file are taken by InfoFasta class
        # You make sure that the arg text is equal to False because the input object is a file and not a list
                
        self.headers = InfoFasta.get_header(self.path_file_isoform, text=False)
        
        
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
            
            
        self.file_random_seq = FileUtils.open_text_a(self.path_random_file)
        
        # The sequences corresponding to header selected are extracted from isoform file   
        for header in random_header:
            Logger.get_instance().info('Header selected : ' + header)
            identifier = header[33:48]
            sequence = InfoFasta.get_seq(self.path_file_isoform, identifier, text=False)
            fasta_seq = SeqToFasta.give_fasta(header, sequence)
            self.file_random_seq.write(fasta_seq)
            
        
        Logger.get_instance().info( " End of selection random sequences \n ")
        
        
    
    # merger_sequences
    # --------------------
    #
    # This method merges the Longest sequences and the random isoform 
    # 
    # The output is the final fasta file of dataset that can be gived as input to DisorderAnalysis
    #
    def merger_sequences(self):
            
        Logger.get_instance().info( " Union of the longest sequences and the random selected isoform ")
            
        # Input variables to merge the longest Novel sequences and random selected isoform of dataset 
        
        self.path_home = Constants.PATH_HOME
        self.path_file_longest = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_FUSION_FILE_LONGEST_PROPERTY, True)
        self.path_file_random = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_FUSION_FILE_RANDOM_PROPERTY, True)
        self.path_final_file = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOWNLOAD_FUSION_FINAL_DATASET_PROPERTY, True)
        
        FileParser.merge_file(self.path_file_longest, self.path_file_random, self.path_final_file)
        
        
        Logger.get_instance().info( " The Final TFs Dataset has been created\n ")


    # whole_procedure
    # ----------------
    # 
    # This method allows to select just some of previously methods
    #
    @staticmethod
    def whole_procedure():
            
        
        # start chrono
        Timer.get_instance().start_chrono()
        Logger.get_instance().info("Start of the sequences extraction.....\n ")
           
        D = DownloadEnsemblSeq()
        
        #D.download_product_gene_seq()
        
        #D.make_dictionary()
        #D.get_longest_seq()
        #D.isoform_sequences()
        #D.merger_sequences()
        
        
            
        Timer.get_instance().stop_chrono(' End of the sequences extraction')
        
        
        
if __name__ == '__main__':
    
    OptionManager.get_instance().initialize()
    
    # Retrieve the  properties DownloadEnsemblSeq
    PropertyManager.get_instance().read_properties( OptionManager.get_instance().get_option( OptionConstants.OPTION_DOWNLOADENSEMBLSEQ_PROPERTY_PATH, True))

    D = DownloadEnsemblSeq()
    
    D.whole_procedure()
    