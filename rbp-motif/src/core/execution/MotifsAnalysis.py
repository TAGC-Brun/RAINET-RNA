# ===================================================================================================
# This class performs the mapping of SLiMs over Iupred,Anchor, DisoRDPbind regions in a given dataset 
# ===================================================================================================


# Starting files and conditions
#    - the iupred, anchor, disoRDPbind analysis have to be already run
#    - the file containing the domain of all proteins in pickle.dump format
#    - file containing the list of proteins
#


from core.util import Constants
from core.util.property.PropertyManager import PropertyManager
from core.data import DataConstants
from core.util.time.Timer import Timer
from core.util.log.Logger import Logger
from core.util.parsing.FileParser import FileParser
from core.util.tools.GlobalOverlapRegionAnalysis import GlobalOverlapRegionAnalysis
from core.util.option import OptionConstants
from core.util.option.OptionManager import OptionManager


class MotifsAnalysis():
    
    def __init__(self):
        
        self.path_home = Constants.PATH_HOME
        self.protein_list_file = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.MOTIFS_PROTEIN_FILE_PROPERTY, True)
        self.protein_list = FileParser.read_file(self.protein_list_file)
        self.motif_folder = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.MOTIFS_FOLDER_PROPERTY, True)
        self.domain_region_file = self.path_home + PropertyManager.get_instance().get_property( DataConstants.MOTIFS_DOMAIN_REGION_FILE_PROPERTY, True)
    
    # iupred_motifs
    # ----------------------------
    #
    # This method calls the GlobalOverlapRegionAnalysis for iupred tool 
    #
    # Output: 
    #    - The result of motifs overlap analysis will be saved in some files in output folder 
    #
    def iupred_motifs(self):
        
        Logger.get_instance().info( "        .....Start of IUPred motifs analysis.....\n")
        
        self.iupred_folder = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.MOTIFS_IUPRED_FOLDER_PROPERTY,True)
        
        # Iupred Analysis at threshold value of 0.4
        
        Timer.get_instance().step(" Start of IUPred motifs analysis - threshold = 0.4 \n")
        
        self.threshold_1 = Constants.MOTIFS_THRESHOLD_1
        self.output_folder_1 = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.MOTIFS_IUP_OUTPUT_FOLDER_1_PROPERTY, True)
        
        GlobalOverlapRegionAnalysis.iupred_overlap_analysis(self.protein_list,self.iupred_folder, self.output_folder_1, self.threshold_1,
                                                            self.motif_folder,self.domain_region_file)
        
        Timer.get_instance().step(" End of IUPred motifs analysis - threshold = 0.4 \n")
        
        
        # Iupred Analysis at threshold value of 0.5
        
        Timer.get_instance().step(" Start of IUPred motifs analysis - threshold = 0.5 \n")
        self.threshold_2 = Constants.MOTIFS_THRESHOLD_2
        self.output_folder_2 = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.MOTIFS_IUP_OUTPUT_FOLDER_2_PROPERTY, True)
        
        GlobalOverlapRegionAnalysis.iupred_overlap_analysis(self.protein_list,self.iupred_folder, self.output_folder_2, self.threshold_2,
                                                            self.motif_folder,self.domain_region_file)
        
        Timer.get_instance().step(" End of IUPred motifs analysis - threshold = 0.5 \n")
        
        Logger.get_instance().info( "        .....End of IUPred motifs analysis.....\n")
    
        
    # anchor_motifs
    # ----------------------------
    #
    # This method calls the GlobalOverlapRegionAnalysis for anchor tool 
    #
    # Output: 
    #    - The result of motifs overlap analysis will be saved in some files in output folder 
    #    
    def anchor_motifs(self):
        
        Logger.get_instance().info( "        .....Start of IUPred motifs analysis.....\n")
        
        self.anchor_folder = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.MOTIFS_ANCHOR_FOLDER_PROPERTY, True)
        self.anchor_output_folder = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.MOTIFS_ANCHOR_OUTPUT_FOLDER_PROPERTY, True)
        
        Timer.get_instance().step(" Start of ANCHOR motifs analysis \n")
        
        GlobalOverlapRegionAnalysis.anchor_overlap_analysis(self.protein_list,self.anchor_folder,self.anchor_output_folder,
                                                            self.motif_folder,self.domain_region_file)
        
         
        Timer.get_instance().step(" End of IUPred motifs analysis \n")
        
        Logger.get_instance().info( "        .....End of ANCHOR motifs analysis.....\n")
        
        
    # disordp_motifs
    # ----------------------------
    #
    # This method calls the GlobalOverlapRegionAnalysis for disordpbind tool 
    #
    # Output: 
    #    - The result of motifs overlap analysis will be saved in some files in output folder 
    
    def disordp_motifs(self):
        
        Logger.get_instance().info( "        .....Start of DisoRDPbind motifs analysis.....\n")
        
        self.disordp_folder = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.MOTIFS_DISORDP_FOLDER_PROPERTY, True)
        self.disordp_output_folder = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.MOTIFS_DISORDP_OUTPUT_FOLDER_PROPERTY, True)
        self.filename =  PropertyManager.get_instance().get_property( DataConstants.MOTIFS_DISORDP_FILE_PROPERTY,True)
        
        Timer.get_instance().step(" Start of DisoRDPbind motifs analysis \n")
        
        GlobalOverlapRegionAnalysis.disordp_overlap_analysis(self.protein_list, self.disordp_folder, self.filename, self.motif_folder,
                                                             self.domain_region_file,self.disordp_output_folder)
        
        
        
        Timer.get_instance().step(" End of DisoRDPbind motifs analysis \n")

        Logger.get_instance().info( "        .....End of DisoRDPbind motifs analysis.....\n")
    
    
    # whole_procedure
    # ----------------
    # 
    # This method allows to select just some of previously methods
    #

    @staticmethod
    def whole_procedure():
            
        dataset_type = PropertyManager.get_instance().get_property( DataConstants.MOTIFS_DATASET_TYPE_PROPERTY, True)
        
        # start chrono
        Timer.get_instance().start_chrono()
        Logger.get_instance().info("        ........Start of " + dataset_type + " Motifs Analysis.....\n ")
        
        motifs = MotifsAnalysis()
        motifs.iupred_motifs()
        motifs.anchor_motifs()
        motifs.disordp_motifs()
        
        Timer.get_instance().stop_chrono(" End of " + dataset_type + " Motifs Analysis")




if __name__ == '__main__':
    
    OptionManager.get_instance().initialize()
    
    # Set the level of verbosity
    Logger.get_instance().set_level( OptionManager.get_instance().get_option( OptionConstants.OPTION_VERBOSITY))

    PropertyManager.get_instance().read_properties( OptionManager.get_instance().get_option( OptionConstants.OPTION_MOTIFS_ANALYSIS_PROPERTY_PATH, True))
    
    MotifsAnalysis.whole_procedure()
    