# ========================================================================
# This class performs the global analysis with regard the Overap between 
# Motifs (SLims) and the region of Iupred or ANCHOR or DisoRDPbind 
# =======================================================================



from core.util.file.FileWriter import FileWriter
from core.util import Constants
from core.util.log.Logger import Logger
from core.util.tools.OverlapMotifsRegions import OverlapMotifsRegions
import math



class GlobalOverlapRegionAnalysis():
    
    #
    # overlap_file_write
    # -------------------
    # 
    # This method allows to perform the OverlapMotifsRegions for a list of proteins 
    # Argument:
    #    - table_tool: containing the slim tool region overlaps and the slim domain overlaps [the table_tool[0] have not to be empty!!!]
    #    - protein: protein analyzed
    #    - toolname: type of toolname
    #    - output folder: folder to save the output data
    # Returns:
    #    - file_slim_region: containing the table of overlap between the slim region and iupred region
    #    - file_domain: containing the table of overlap between the slim region and iupred region
    #                    only if the slim fall into a domain
    @staticmethod
    def overlap_file_write(table_tool,output_folder,protein,toolname):
        
        dict_filename = Constants.DICT_FILE_OUTPUT
        
        filename_slim_region =  dict_filename[toolname][0] + protein + Constants.EXTENSION_TXT_FILE
        filename_doamin_region = dict_filename[toolname][1] + protein + Constants.EXTENSION_TXT_FILE
        filepath_slim = output_folder + filename_slim_region
        filepath_domain = output_folder + filename_doamin_region
        
        title_slim_region = [[Constants.SLIM_NAME,Constants.SLIM_REGION_COUNTER,Constants.START_SLIM_REGION,Constants.END_SLIM_REGION,Constants.TOOLNAME,
                              Constants.TOOLNAME_REGION_COUNTER,Constants.START_TOOL_REGION, Constants.END_TOOL_REGION, 
                              Constants.OVERLAP_OUTCOME, Constants.OVERLAP_LENGTH]] 
        
        table_slim_region = title_slim_region + table_tool[0]
        table_domain_region = table_tool[1] 

        if table_tool[0] != []:
            FileWriter.write_table(filepath_slim, table_slim_region)
        else:
            Logger.get_instance().debug("  The slim overlap file " + filename_slim_region+ " has not been written")
        if table_tool[1] != []:
            FileWriter.write_table(filepath_domain, table_domain_region)
        else:
            Logger.get_instance().debug("  The domain overlap file " + filename_doamin_region+ " has not been written\n\n")


    #
    # iupred_overlap_analysis
    # -------------------
    # 
    # This method allows to perform the OverlapMotifsRegions for a list of proteins 
    # Argument:
    #    - protein_list: list of protein
    #    - iupred folder: directory containing the tool prediction
    #    - motif_folder: directory containing the SLiM Found with ANCHOR tool
    #    - output folder: folder to save the output data
    #    - threshold: threshold to select the iupred region
    # Returns:
    #    the following file for all proteins in protin_list
    #    - file_slim_region: containing the table of overlap between the slim region and iupred region
    #    - file_domain: containing the table of overlap between the slim region and iupred region
    #                    only if the slim fall into a domain
    @staticmethod
    def iupred_overlap_analysis(protein_list,iupred_folder, output_folder, threshold, motif_folder,domain_region_file):
        
        toolname = Constants.IUPRED_TOOL

        for num, protein in enumerate(protein_list):
            prot_counter = num +1
            result = math.fmod(prot_counter,100)
            if result == 0.0:
                Logger.get_instance().info('        ' +str(prot_counter)+ "  Protein analyzed: " + protein) 
            # This line takes the iupred region of a protein
            iupred_region =  OverlapMotifsRegions.iupred_region(iupred_folder,threshold, protein )
            # This line calls the method in order to check the overlap 
            table_iupred = OverlapMotifsRegions.overlap_analysis(protein,iupred_region,toolname, motif_folder,domain_region_file)
            # table Writing
            GlobalOverlapRegionAnalysis.overlap_file_write(table_iupred, output_folder, protein, toolname)

    
    #
    # anchor_overlap_analysis
    # -------------------
    # 
    # This method allows to perform the OverlapMotifsRegions for a list of proteins 
    # Argument:
    #    - protein_list: list of protein
    #    - anchor folder: directory containing the tool prediction
    #    - motif_folder: directory containing the SLiM Found with ANCHOR tool
    #    - output folder: folder to save the output data
    # Returns:
    #    the following file for all proteins in protin_list
    #    - file_slim_region: containing the table of overlap between the slim region and iupred region
    #    - file_domain: containing the table of overlap between the slim region and iupred region
    #                    only if the slim fall into a domain
    @staticmethod
    def anchor_overlap_analysis(protein_list,anchor_folder,output_folder,motif_folder,domain_region_file):
        
                
        toolname = Constants.ANCHOR_TOOL
        
        for num, protein in enumerate(protein_list):
            prot_counter = num +1
            result = math.fmod(prot_counter,100)
            if result == 0.0:
                Logger.get_instance().info('        ' + str(prot_counter)+ "th  Protein analyzed: " + protein) 
            # This line takes the iupred region of a protein
            anchor_region = OverlapMotifsRegions.anchor_region(anchor_folder, protein )
            table_anchor = OverlapMotifsRegions.overlap_analysis(protein,anchor_region,toolname, motif_folder,domain_region_file )
            # table Writing
            GlobalOverlapRegionAnalysis.overlap_file_write(table_anchor, output_folder, protein, toolname)

    #
    # disordp_overlap_analysis
    # -------------------
    # 
    # This method allows to perform the Overlapanalysis for a list of proteins 
    # Argument:
    #    - protein_list: list of protein
    #    - disordp folder: directory containing the tool prediction file
    #    - filename: file containing the disordpb prediction
    #    - motif_folder: directory containing the SLiM Found with ANCHOR tool
    #    - output folder: folder to save the output data
    # Returns:
    #    the following file for all proteins in protin_list
    #    - file_slim_region: containing the table of overlap between the slim region and iupred region
    #    - file_domain: containing the table of overlap between the slim region and iupred region
    #                    only if the slim fall into a domain
    @staticmethod
    def disordp_overlap_analysis(protein_list, disordp_folder, filename, motif_folder,domain_region_file,output_folder):
        
        
        toolname = Constants.DISORDP_TOOL
        
        
        for num, protein in enumerate(protein_list):
            prot_counter = num +1
            result = math.fmod(prot_counter,100)
            if result == 0.0:
                Logger.get_instance().info('        ' +str(prot_counter)+ " Protein analyzed: " + protein)
            # This line takes the iupred region of a protein
            disordp_region = OverlapMotifsRegions.disordp_region(disordp_folder, filename, protein )
            table_disordp = OverlapMotifsRegions.overlap_analysis(protein,disordp_region, toolname, motif_folder,domain_region_file )
            # table Writing
            GlobalOverlapRegionAnalysis.overlap_file_write(table_disordp, output_folder, protein, toolname)

        
        
        
        