#===============================================================================
#===============================================================================
# COMPARISON BETWEEN TWO DATASET
#===============================================================================
#===============================================================================


from core.util.parsing.FileParser import FileParser
from core.util.parsing.TableWrapper import TableWrapper
from core.util.file.FileWriter import FileWriter
from core.util import Constants
from core.util.log.Logger import Logger
from core.util.option.OptionManager import OptionManager
from core.util.option import OptionConstants





class InfoDataset():
    
    # get_dataset_feature
    # ----------------------
    #
    # This method extracts one features of a datasets and make a list
    # 
    # Arguments:
    #    - dataset_path: path of dataset file
    #    - col_feature: number of column feature
    #    - skip: number of line to skip. Default value is True (it is presumed that the dataset have a header line)
    # Returns:
    #    - feature that you want extract from dataset
    #
    @staticmethod
    def get_dataset_feature(dataset_path, col_feature, skip=1):
        dataset = FileParser.make_table(dataset_path, skip=skip)
        feature = TableWrapper.get_column(dataset, col_feature)
        return feature
    
    # check_length
    # ----------------------
    #
    # This method checks if the length of feature is that you expected
    #
    # Arguments:
    #    - feature: general column of which you want check the length
    #    - expected_length
    # Returns:
    #    length: length of selected feature 
    #
    @staticmethod
    def check_length(feature, expected_length=None):

        length = len(feature)
        if expected_length is not None:
            if length == expected_length:
                Logger.get_instance().info(" The length of dataset feature is that you expected\n")
            else:
                Logger.get_instance().info(" The length of dataset feature isn't that you expected\n")
        else:
            pass
        return length
    
    
    # comparison_dataset
    # ------------------------------
    # 
    # This method compares a feature of two dataset and it returns a tuple 
    # containing (intersection, difference_a-b, difference_b-a, union)
    #
    # Arguments:
    #    -f_dataset_a: feature dataset 1
    #    -f_dataset_b: feature dataset 2
    # Returns:
    #    -tuple: common item, diff_ab item, diff_ba item and union of items
    #
    @staticmethod
    def comparison_dataset(f_dataset_a, f_dataset_b):


        a = set(f_dataset_a)
        b = set(f_dataset_b)
        common = a & b
        diff_ab = a - b
        diff_ba = b - a
        union = a | b
        return (common, diff_ab, diff_ba, union)
    
    
    # global_analysis_dataset
    # -----------------------
    #
    # This method makes a global analysis starting from two datasets
    # 
    # Arguments:
    #    - path_input: datasets directory 
    #    - dataset: tuple object with file name (dataset_a, dataset_b)
    #    - index_col: tuple containing the i-th column index of datasets
    #    - path_output: directory path of output file
    #    - length: tuple containing expected_len_a, expected_len_b
    #
    # Returns:
    #    - files containing the common and different item of datasets
    #    - Note: the different item are dataset b - dataset a
    #
    @staticmethod
    def global_analysis_dataset(path_input, dataset, index_col, path_output, length=None ):
        
        (dataset_a, dataset_b) = dataset
        col_a = int(index_col[0])
        col_b = int(index_col[1])
        
        path_dataset_a = path_input + dataset_a
        path_dataset_b = path_input + dataset_b
        
        # dataset a
        # ==========================================================
        feature_a = InfoDataset.get_dataset_feature(path_dataset_a, col_a)
        
        length_a = InfoDataset.check_length(feature_a)
        Logger.get_instance().info(' The length of dataset a is : ' + str(length_a))
                                   
        # dataset 2  
        # ==========================================================
        feature_b = InfoDataset.get_dataset_feature(path_dataset_b, col_b)
        
        length_b = InfoDataset.check_length(feature_b)
        Logger.get_instance().info(' The length of dataset b is : ' + str(length_b))
        
        
        if length is None:
            pass
        else:
            expected_len_a = int(length[0])
            expected_len_b = int(length[1])
            Logger.get_instance().info('Dataset 1')
            length_a = InfoDataset.check_length(feature_a, expected_length=expected_len_a)
            Logger.get_instance().info('Dataset 2')
            length_b = InfoDataset.check_length(feature_b, expected_length=expected_len_b)   

            
        
        
        operations = InfoDataset.comparison_dataset(feature_a, feature_b)
        common_items = operations[0]
        diff_item_ab = operations[1]
        diff_item_ba = operations[2]
        
        # Common items
        Logger.get_instance().info(' The number of common items to a and b dataset :\
' + str(len(common_items)))
        
        # Items in b and not in a
        Logger.get_instance().info(' The number of items that are present \
in dataset b but not in dataset a :' + str(len(diff_item_ba)))
        
        # Items in a and not in b
        Logger.get_instance().info(' The number of items that are present \
in dataset a but not in dataset b :' + str(len(diff_item_ab)))
        
        
        # Writing dataset intersection and difference
        # ==========================================================
        
        namefile_common = Constants.FILE_COMMON
        common_file = path_output + namefile_common 
        FileWriter.write_table(common_file, common_items)

        namefile_difference = Constants.FILE_DIFF
        diff_file = path_output + namefile_difference
        FileWriter.write_table(diff_file, diff_item_ba)
        
        
        

if __name__ == '__main__':
    
    OptionManager.get_instance().initialize()
    
    InfoDataset.global_anlysis_dataset(OptionManager.get_instance().get_option( OptionConstants.OPTION_INPUT_PATH),
                                       OptionManager.get_instance().get_option( OptionConstants.OPTION_MORE_FILENAME),
                                       OptionManager.get_instance().get_option( OptionConstants.OPTION_NUMBER_INDEX),
                                       OptionManager.get_instance().get_option( OptionConstants.OPTION_OUTPUT_PATH),
                                       length=OptionManager.get_instance().get_option( OptionConstants.OPTION_NUMBER_LENGTH),
                                       )
    
    

