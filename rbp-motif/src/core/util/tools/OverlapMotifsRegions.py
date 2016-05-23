





from core.util.parsing.FileParser import FileParser
from core.util import Constants
from core.util.tools.Anchor import Anchor
from core.util.tools.Iupred import Iupred
from core.util.tools.DisoRDPbind import DisoRDPbind
from core.util.log.Logger import Logger
from core.util.tools.OverlapTwoRegion import OverlapTwoRegion
from core.util.tools.OverlapSlimDomain import OverlapSlimDomain
import os.path


class OverlapMotifsRegions():
    
    
    # elm_region
    # -------------------
    # 
    # Argument:
    #    - protein name
    # Returns:
    #    - dictionary containing as key the motifs names and as value the regions[start end]
    
    @staticmethod
    def elm_region(protname, directory):
        

        filename = Constants.ELM_FILE + protname + Constants.EXTENSION_TXT_FILE
        filepath = directory + filename
        dictionary = {}
        
        occurrence = os.path.isfile(filepath)

        if occurrence == True:
            motif_table = FileParser.make_table(filepath, skip=1)
            # for loop to make the dictionary
            for m, row in enumerate(motif_table):
                name = row[0]
                start = int(row[1])
                end = int(row[2])
                if name not in dictionary:
                    dictionary[name] = [[start, end]]
                else:
                    dictionary[name].append([start,end])
            Logger.get_instance().info(" In the protein " + protname +' '+ str(len(dictionary)) + ' SLiMs have been found')      
        else:
            Logger.get_instance().info(' This protein has not SLiMs ')
        return dictionary
            
            

        
    
    
    # iupred_region
    # -------------------
    # 
    # Argument:
    #    - protein name: ENSEMBL identifier
    #    - input_folder: directory containing the Iupred output
    #   
    # Returns:
    #    - iupred_region_table : table containing the disordered region table

        
    @staticmethod    
    def iupred_region(iupred_folder,threshold, protname ):
        
        
        filename = Constants.IUPRED_PRED_FILE + protname + Constants.EXTENSION_TXT_FILE

        filepath = iupred_folder + filename 

        # extraction of Iupred regions 
        # the information are memorized in a table
        dict_info = Iupred.iupred_info(filepath, protname, float(threshold))
        # Table containing the disordered regions 
        key_dictionary = Constants.KEY_IUPRED
        iupred_table = dict_info[key_dictionary]
        iupred_region_table = [ [int(line[2]), int(line[3])] for line in iupred_table]
    
        return iupred_region_table
    
    # anchor_region
    # -------------------
    # 
    # Argument:
    #    - protein name: ENSEMBL identifier
    #    - anchor_folder: directory containing the ANCHOR output
    # 
    # Returns:
    #    - anchor_region_table : table containing the binding region table


        
    @staticmethod    
    def anchor_region(anchor_folder, protname):
        
        filename = Constants.ANCHOR_PRED_FILE + protname + Constants.EXTENSION_TXT_FILE

        filepath = anchor_folder + filename 

        # extraction of anchor regions 
        # the information are memorized in a table
        dict_info = Anchor.anchor_info(filepath, protname)
        key_dictionary = Constants.KEY_ANCHOR
        # Table containing the anchor regions 
        anchor_table = dict_info[key_dictionary]
        anchor_region_table = [ [int(line[2]), int(line[3])] for line in anchor_table]
    
        return anchor_region_table      
         
    
        
    # disordbp_region
    # -------------------
    # 
    # Argument:
    #    - protein name: ENSEMBL identifier
    #    - directory: folder containing the disordp_file
    #    - disordp_file: name of file containing the output of disoRDPbind tool
    # 
    # Returns:
    #    - disordp_region_table : table containing the RNA-binding region table


        
    @staticmethod    
    def disordp_region(directory, disordp_file, protname):
        

        filepath = directory + disordp_file 
        
        # the information about disordbp bind are memorized in a file containing information of many proteins
        # this command allows to select the output information about one only protein (protname)
        output_proteins = DisoRDPbind.output_reading(filepath)
        
        proteins = [ line.split('\n')[0] for line in output_proteins]
        if protname in proteins:
            prot_selected = [ line for line in output_proteins if line.split('\n')[0]==protname ]
            if 'WARNING:' in prot_selected[0]:
                Logger.get_instance().warning( " This protein contains >=10000 residues\
 (DisoRBDbind cannot predict the proteins with size >=10000) " + protname)
                disordp_region_table = []
            else:   
                # extraction of disordp regions 
                # the information are memorized in a table
                dict_info = DisoRDPbind.fraction_calculation(prot_selected[0])
                key_dictionary = Constants.KEY_DISORDP
                # Table containing the anchor regions 
                disordp_table = dict_info[key_dictionary]
                disordp_region_table = [ [int(line[2]), int(line[3])] for line in disordp_table]
            
        else:
            Logger.get_instance().warning(' This protein is not in DisoRDPbind prediction \
(DisoRBDbind cannot predict the proteins with size of 4 amino acids) ' + protname)
            disordp_region_table = []
            
        return disordp_region_table      
         
    
    
    
    # overlap_anlysis
    # -------------------
    #
    # Arguments
    #     - motif_directory: folder containing the motif found
    #     - domain_region_file: dictionary (pickle dump) file containing the 
    #       domains of proteins
    # Returns
    #     - tuple: containing two table 
    #            table 1 contains the comparison between slim and disordered regions
    #            table 2 contains the comparison between slim and disordered regions 
    #            if the slim occurs in a domain regions  
    #            
    
    
    @staticmethod    
    def overlap_analysis(protname,tool_table,toolname, motif_directory,domain_region_file ):
        
        Logger.get_instance().info(' The protein analyzed is '+ protname)
        
        # SLiMs regions 
        # dictionary containing the slims and regions as value for each key
        elm_region = OverlapMotifsRegions.elm_region(protname, motif_directory)
        
        # Protein Domain
        prot_domain = OverlapSlimDomain.domain_one_protein(domain_region_file, protname)
        num_domain = len(prot_domain)
        Logger.get_instance().info(' Protein domains: '  + str(num_domain))

        
        
        # Tool regions 
        num_tool_region = len(tool_table)
        Logger.get_instance().info(' '+ toolname + ' regions: ' +str(num_tool_region) +'\n' )
        
        if len(elm_region.keys()) == 0 or num_tool_region == 0:
            return ([],[])
        else:            
            # This section aims to find if there are Motifs in some specific regions
            
            table_all_slim_region = []
            table_overlap_domain = []
            table_domain = []
            # the first for loop flows over the slim contained in elm_region 
            # the second for loop takes one region of slim at time and check the possible
            # overlap with the domain of proteins and the third for loop 
            # flows over the region of a tool(iupred,anchor,disoRDPbind) in order to see if there are overlap between the 
            # regions and the slim motif region
            # ---------------------------------
            # k_slim counter: counts the number of slim in the protein
            # r_slim counter: counts the number of slim regions 
            # r_tool counter: counts the number of tool regions
            for k_slim, slim in enumerate(sorted(elm_region.keys())):
                total_slim = len(elm_region.keys())
                Logger.get_instance().debug(' The SLiM analyzed is ' + slim +' (' + str(k_slim+1)+' out of ' + str(total_slim)+ ')' )
                Logger.get_instance().debug(' The total occurrences number of this SLiMs is  '+ str(len(elm_region[slim])) )
                Logger.get_instance().debug(" ----------------------------------------------- \n")
                for r_slim, slim_region in enumerate(elm_region[slim]):
                    table_slim_region = []
                    if num_domain != 0:
                        Logger.get_instance().debug('    - '+  str(r_slim+1) +'th SLiM region analyzed ('+ protname+ ')\n') 
                        # Comparison between slim region and domain region
                        table_domain = OverlapSlimDomain.overlap_motif_domain(domain_region_file, protname, slim, slim_region)
                    else:
                        pass
                    # Comparison between Slim region and tool region
                    Logger.get_instance().debug(' Comparison of SLiM region and ' + toolname+ ' regions ')
                    for r_tool, region in enumerate(tool_table):
                        Logger.get_instance().debug('            - '+  str(r_tool+1) +'th '+toolname +' region analyzed') 
                        # Calling of function that check if two generic region are overlapped
                        # return the position of overlap into list if the regions are not overlapped
                        # the list is empty
                        overlap_region = OverlapTwoRegion.region_overlap(slim_region, region)
                        # Definition of start and end region position
                        start_elm = slim_region[0]
                        end_elm = slim_region[1]
                        start_iupred = region[0]
                        end_iupred = region[1]
                        length_overlap = len(overlap_region)
                        # if the regions are overlapped the overlap region variable will be long at least one item
                        if len(overlap_region) >= 1:
                            overlap = Constants.OVERLAP_YES
                        else:
                            overlap = Constants.OVERLAP_NO
                        # Definition of row that is saved in a table in order to be memorized in a file successively
                        row_slim = [slim, str(r_slim+1), start_elm, end_elm, toolname, r_tool+1, start_iupred, end_iupred, overlap, length_overlap]
                        Logger.get_instance().debug("        Overlap result between " + slim+ " and " + str(r_tool+1) +"th region: " ' Overlap? '+ overlap+' \n' )
                        # table containing the overlap between one slim and all tool regions
                        table_slim_region.append(row_slim)
                    Logger.get_instance().debug("                 /-----------------/\n ") 
                    # If there is an overlap between one motifs and a domain the information are saved in
                    # a different table in order to keep all information but separated
                    if len(table_domain) >= 1:
                        table_overlap_domain = table_overlap_domain + table_domain + table_slim_region
                    # when the overlap between slim and motifs is not occurred the main table is expanded 
                    # with new overlap information
                    else:
                        table_all_slim_region = table_all_slim_region  + table_slim_region
                        
            # the output is composed of two table that will be used to write a file            
            return (table_all_slim_region, table_overlap_domain)
    
    
    
    
               
        
        
        