# ==============================================================
# This class check if two region are overlapped or not
# =============================================================


import pickle
from core.util.file.FileUtils import FileUtils
from core.util.tools.OverlapTwoRegion import OverlapTwoRegion
from core.util import Constants
from core.util.log.Logger import Logger
from core.util.parsing.TableWrapper import TableWrapper

class OverlapSlimDomain():
    
    # domain_all_protein
    # -------------------
    # 
    # Arguments:
    #     - domain_region_file : dictionary containing the protein as key and its 
    #      own domain as value in dictionary format the domains are the keys and 
    #      the values the regions
    # Returns:
    #     - dictionary uploaded
    

    
    @staticmethod
    def domain_all_protein(domain_region_file):
        
        file_domain = FileUtils.open_text_r(domain_region_file)
        
        # Importation of Dictionary 
        dict_domain = pickle.load(file_domain)
        
        return dict_domain
    
    # domain_one_protein
    # -------------------
    # 
    # Arguments:
    #     - domain_region_file : dictionary containing the protein as key and its 
    #      own domain as value in dictionary format the domains are the keys and 
    #      the values the regions
    #     - protname
    # Returns:
    #     - dictionary containing the domain of specific proteins
    
    @staticmethod
    def domain_one_protein(domain_region_file, protname):
        
        file_domain = FileUtils.open_text_r(domain_region_file)
        
        # Importation of Dictionary 
        dict_domain = pickle.load(file_domain)
        
        if protname in dict_domain:
            domain_prot = dict_domain[protname]
            return domain_prot
        else:
            Logger.get_instance().debug(" Protein without domains " + protname)
            return []

        
    
    # regios_overlap
    # -------------------
    # 
    # Comparison of the possible overlap between slim and domain
    # Argument:
    #   
    #    - protname: identifier of protein
    #    - slim_name: slim that you want compare
    #    - slim_region: region of slim   
    # Returns:
    #    - table: containing the information of overlap between regions if there not are
    #             overlap the table is empty
    #
    
    @staticmethod
    def overlap_motif_domain(domain_region_file, protname, slim_name, slim_region):
        
        # uploading of dictionary containing of protein domains
        prot_domain = OverlapSlimDomain.domain_one_protein(domain_region_file, protname)
        
        # Selection of proteins domain
        if prot_domain != {}:
            # for loop to find the possible overlap between domain and slim
            table_domain = []
            Logger.get_instance().debug(" Regions comparison between " +slim_name+' and protein domains') 
            Logger.get_instance().debug(" ----------------------------------------------------------------- \n")
            for n_domain, domain in enumerate(sorted(prot_domain)):
                domain_regions = prot_domain[domain] # table
                Logger.get_instance().debug( '        - '+ str(n_domain+1) +'th domain analyzed '+ domain)
                for n_reg, region_dom in enumerate(domain_regions):
                    Logger.get_instance().debug( '        - '+ str(n_reg+1) +'th region analyzed of domain: '+domain)
                    overlap_slim_domain = OverlapTwoRegion.region_overlap(slim_region, region_dom)
                    if len(overlap_slim_domain) >= 1:
                        overlap_mot_domain = Constants.OVERLAP_YES
                        row_domain = [slim_name, 0, slim_region[0], slim_region[1], domain, 0, region_dom[0],region_dom[1],overlap_mot_domain, len(overlap_slim_domain)]
                        table_domain.append(row_domain)
                    else: 
                        overlap_mot_domain = Constants.OVERLAP_NO
                Logger.get_instance().debug("        Overlap result between " + slim_name+ " and " + domain +": " + ' Overlap? '+ overlap_mot_domain+' \n' )    
            
            num_domain = len(table_domain)
            if num_domain > 1:
                domain = TableWrapper.get_column(table_domain, 4)
                region_start =  TableWrapper.get_column(table_domain, 6)
                region_end = TableWrapper.get_column(table_domain, 7)
                region = [ (str(start), str(end)) for start, end in zip(*[region_start,region_end])]
                new_regions = ['-'.join(reg) for reg in region]
                string_slim_region = [str(item) for item in slim_region]
                Logger.get_instance().debug(" The motif: "+slim_name+ ' ('+'-'.join(string_slim_region)+") falls into at least two domain instances or two different domains: " 
                                           + ','.join(domain)+' '+ ','.join(new_regions)+'\n')
            return table_domain
            
        else:
            return []