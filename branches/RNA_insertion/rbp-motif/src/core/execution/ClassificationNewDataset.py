# =====================================================================================
# This class has been created in order to classify the RBP dataset in according to:
# A) RNA PUTATIVE TARGET
# B) CLASS DOMAIN
# =====================================================================================


# Starting informations
#
# A) RNA PUTATIVE TARGET
#
# The RBP dataset is not homogeneous with regard the information about RNA Putative Target
# This is due to the fact that the RBP dataset is made starting from two different datasets
# 
# Dataset 1: each protein has the corresponding RNA putative target
# Dataset 2: the proteins haven't the corresponding RNA putative target but this information
#            can be deduced considering other informations provided by dataset
#            Indeed the dataset provides the information about the protein attendance in other
#            datasets of which is known the putative RNA target 
# Other datasets:
#     castello --> mRNA
#     baltz --> mRNA
#     esc --> mRNA
#     rbpdb --> RBP protein --manually check of putative RNA targrt
#     rnacompete --> Unkonwn
#
# B) CLASS DOMAIN
#
# The second dataset provides a list of domain that are classified in:
#     - Classical
#     - Non-Classical
#     - Unclassified (Unknown)
# 
# The purpose is to classify the protein according to type of domain in the whole RBP dataset
# It need to take in account that the datasets have a different nomenclature domains that reflects on 
# file parsing. 
# 
#


from core.util import Constants
from core.util.property.PropertyManager import PropertyManager
from core.data import DataConstants
from core.util.format.InfoFasta import InfoFasta
from core.util.file.FileUtils import FileUtils
from core.util.parsing.FileParser import FileParser
from core.util.parsing.TableWrapper import TableWrapper
from core.util.file.FileWriter import FileWriter



class ClassificationNewDataset():
    
    # A) RNA PUTATIVE TARGET
    
    # mrna_protein
    # ------------------------
    #
    # This method allow to identify the protein of dataset 2 (only different gene) that have as putative RNA target the mRNA
    #     
    def mrna_protein(self): 
    
        self.path_home = Constants.PATH_HOME
        self.jproteomics_seq =  self.path_home + PropertyManager.get_instance().get_property( DataConstants.PUTATIVERNA_JPROTEOMICS_SEQ_PROPERTY, True)
        self.jproteomics_info = self.path_home + PropertyManager.get_instance().get_property( DataConstants.PUTATIVERNA_JPROTEOMICS_INFO_PROPERTY, True)
        
        # Reading J proteomics file containing information about the protein attendance in others datasets
        #
        info_table = FileParser.make_table(self.jproteomics_info, '\t', skip=1)
        #
        # column Extraction
        gene_info = TableWrapper.get_column(info_table, 1)
        castello = TableWrapper.get_column(info_table, 2)
        baltz = TableWrapper.get_column(info_table, 3)
        esc = TableWrapper.get_column(info_table, 4)
        rbpdb = TableWrapper.get_column(info_table, 5)
        rnacompete = TableWrapper.get_column(info_table, 6)
        
        mrna = [(cast_val, baltz_val, esc_val) for cast_val, baltz_val, esc_val in zip(castello,baltz,esc)]
        
        self.header = InfoFasta.get_header(self.jproteomics_seq, text=False)
         
        # Protein and gene of dataset 2 attend in final RBP dataset
        # 
        gene_seq = [ item.split('>')[1].split('|')[0] for item in self.header]
        prot_seq = [ item.split('>')[1].split('|')[2] for item in self.header]
        
        
        # Construction of column containing the putative RNA target for each gene 
        # in according to if condition    
        putative_rna_target = []
        rbpdb_target = [] 
        for n, gene_id in enumerate(gene_seq):
            print gene_id
            ind_gene_id = gene_info.index(gene_id)
            # if the gene have at least one Y in this array can be considered to have mRNA target
            if 'Y' in mrna[ind_gene_id]:
                rna_target = 'mRNA'
                print 'The putative Rna target of this gene is ' +  rna_target
                putative_rna_target.append([gene_id, prot_seq[n], rna_target])
            elif 'N' in mrna[ind_gene_id]:
                # if the gene attend in rnacompete means that the RNA target is unkown
                if 'N' in rbpdb[ind_gene_id] and 'Y' in rnacompete[ind_gene_id]:
                    rna_target = 'unknown'
                    print 'The putative Rna target of this gene is ' +  rna_target
                    putative_rna_target.append([gene_id, prot_seq[n], rna_target])
                # if the gene attend in rbpdb means just it is a RBP protein (no information about RNA target)
                elif 'Y' in rbpdb[ind_gene_id] and 'N' in rnacompete[ind_gene_id]:
                    rna_target = 'RBPDB'
                    print 'The gene is in ' +  rna_target
                    rbpdb_target.append([gene_id, prot_seq[n], rna_target])
            else:
                print ' Check this line' + info_table[ind_gene_id]
                
        
        
        self.file_rna_target_jeproteomics = self.path_home + PropertyManager.get_instance().get_property( DataConstants.PUTATIVE_MRNA_GENE_JPROTEOMICS_PROPERTY, True)
        self.file_rbpdb_jproteomics = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.PUTATIVE_RBPDB_GENE_JPROTEOMICS_PROPERTY, True)
        
        # File Writing
        
        FileWriter.write_table(self.file_mrna_jeproteomics, putative_rna_target , symbol='\t')
        
        FileWriter.write_table(self.file_rbpdb_jproteomics, rbpdb_target , symbol='\t')                                                               
    
                  
    # rna_target
    # ------------------
    #               
    # This method divided the RBP of first dataset (NatRevGentics) in according to type of putative RNA target
    # 
    # This method produces as many proteins file as the type of putative RNA are
    def rna_target(self):
        
        
        self.path_home = Constants.PATH_HOME
        self.file_seq_natrevgenetics = self.path_home + PropertyManager.get_instance().get_property( DataConstants.PUTATIVERNA_NATREVGENETICS_SEQ_PROPERTY, True)
        self.natrevgenetics_info = self.path_home + PropertyManager.get_instance().get_property( DataConstants.PUTATIVE_RNA_NATREVGENETICS_INFO_PROPERTY, True)
        
        
        info_table = FileParser.make_table(self.natrevgenetics_info, '\t', skip=1)
        
        prot_info = TableWrapper.get_column(info_table, 1)
        putative_rna = TableWrapper.get_column(info_table, 3)
        
        
        self.header = InfoFasta.get_header(self.file_seq_natrevgenetics, text=False)
        
        gene_seq = [ item.split('>')[1].split('|')[0] for item in self.header]
        prot_seq = [ item.split('>')[1].split('|')[2] for item in self.header]
        
    
        # Creation of Table containing gene id, prot id and rna target
        
        putative_rna_target = []
        type_rna_target = []
        for n, prot in enumerate(prot_seq):
            print prot
            index_prot = prot_info.index(prot)
            rna_target = putative_rna[index_prot]
            row = [gene_seq[n], prot_seq[n], rna_target]
            type_rna_target.append(rna_target)
            putative_rna_target.append(row)
            print " The putative rna target is " + rna_target
            

            
        self.file_all_rna_target = self.path_home + PropertyManager.get_instance().get_property( DataConstants.PUTATIVE_ALL_RNA_TARGET_PROPERTY, True)
        
        
        FileWriter.write_table(self.file_all_rna_target, putative_rna_target , symbol='\t')
        
        # set of RNA target type in order to create different list 
        # 
        unique_rna_target = set(type_rna_target)
        
        
        info_new_table = FileParser.make_table(self.file_all_rna_target, '\t')
        
        # Columns extraction
        prot_name = TableWrapper.get_column(info_new_table, 1)
        type_rnatarget = TableWrapper.get_column(info_new_table, 2)
        
        file_output = self.path_home +  PropertyManager.get_instance().get_property( DataConstants.PUTATIVE_RNA_OUTPUT_PROPERTY, True)
        
        # this for loop allows to create a proteins files for each RNA target type 
        for item in unique_rna_target:
            file_name = file_output + item + PropertyManager.get_instance().get_property( DataConstants.PUTATIVE_RNA_TARGET_DATASET_NAME_PROPERTY, True)
            file_rna = FileUtils.open_text_a(file_name)
            for n, type_rna in enumerate(type_rnatarget):
                if type_rna == item:
                    file_rna.write(prot_name[n])
            
            file_rna.close()
                
    
    
    # B) DOMAIN CLASSIFICATION
    
    
    # domain_classification
    # --------------------------------
    #      
    # The domain classification is performed for all protein of dataset 1 and dataset 2
    # At the end a table is created this contains the gene name, protein name, domain, pfam id, class domain
    #    
    def domain_classification(self):
        
        self.path_home = Constants.PATH_HOME
        self.file_domain = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOMAIN_LIST_FILE_PROPERTY, True)
        self.file_jprot_information = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOMAIN_LIST_JPROTEOMICS_PROPERTY, True)
        self.file_pfamid = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOMAIN_FILE_PFAM_PROPERTY, True)
        
        # reading of pfam table in particular of the pfam id and the domain name                                                                                 
        table_pfam = FileParser.make_table(self.file_pfamid)
        pfamid = TableWrapper.get_column(table_pfam, 0)
        domain_name = TableWrapper.get_column(table_pfam, 3)
        
        # dictionary pfamid and domain name        
        dict_pfam = { row[3]: row[0] for row in table_pfam }
        dict_pfam['DUF1785'] = '-'
        dict_pfam['DUF1898'] = '-'
        
        # reading of domain classification provided by dataset 2
        # dictionary motifs--> class
        
        table_domain = FileParser.make_table(self.file_domain)
        dict_type_domain = TableWrapper.make_dictionary(table_domain)
                
        # make a inverse dictionary class--> motifs
        inverse_table_domain = TableWrapper.inv_column(table_domain, 0, 1)
        dict_class_domain = TableWrapper.make_dictionary(inverse_table_domain)
        
        
        # reading of domain information of dataset 2
        jprot_table = FileParser.make_table(self.file_jprot_information, skip=1)
        
        for i, item in enumerate(jprot_table):
            if len(item) == 1:
                jprot_table[i] = [item[0], '.', '.']
                
                
                
        # Jproteomics
        # ===================================

        # extraction of domain for each gene
        # extraction of gene id
        domain_column_jprot = TableWrapper.get_column(jprot_table, 1)
        genes_jprot = TableWrapper.get_column(jprot_table, 0)
        pfam_id_jprot = TableWrapper.get_column(jprot_table, 2)
       

        
        # This part reads the type of domain for each gene and creates a new column with the type of domain
        # at the end returns a new table with namegene, domain name, pfam id, class of domain
        #
        # One gene/protein can have more than one domain
        # For each gene the following part checks if the classification of protein domains by comparing the domain with dict_class_domains:
        # 
        # 
        
        count_classical = 0
        count_nonclassic = 0
        count_unclissified = 0
        count_any_class = 0
        count_no_domain = 0
        new_table_jprot = []
        for i, gene in enumerate(genes_jprot):
            row = []
            row.append(gene)
            row.append(domain_column_jprot[i])
            row.append(pfam_id_jprot[i])
            string_domains = ''
            # If the protein hasn't any domain, add static 'no-domains' information
            if domain_column_jprot[i] == '.':
                string_domains += 'No-domains,'
                count_no_domain += 1
            # If the protein contains some domains checks the class and 
            # makes a string containing the class domains separated by a comma
            else: 
                domains = domain_column_jprot[i].split(',')
                for type_domain in domains:
                    print type_domain
                    if type_domain in dict_type_domain:
                        class_domain = dict_type_domain[type_domain][0]
                        if class_domain == 'Classical':
                            string_domains += 'Classical,'
                            count_classical += 1
                        elif class_domain == 'Non-classical':
                            string_domains += 'Non-classical,'
                            count_nonclassic += 1
                        elif class_domain == 'Unclassified':
                            string_domains += 'Unclassified,'
                            count_unclissified += 1
                    elif type_domain not in dict_type_domain:
                        string_domains += 'Other-Domain,'
                        count_any_class += 1                        
                    else:
                        print 'unexpected case', type_domain
            # -1 allows to delete the last comma in string
            row.append(string_domains[0:len(string_domains)-1])
            new_table_jprot.append(row)
                    
        
        # print of proteins number for each domain class
        print count_classical
        print count_nonclassic
        print count_unclissified
        print count_any_class
        print count_no_domain
            
        self.path_ouput_file_jprot = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOMAIN_FINAL_TABLE_JPROT_PROPERTY, True)
        FileWriter.write_table(self.path_ouput_file_jprot, new_table_jprot)
            
        
        
        
        # NatRevGenetics
        # =========================================================
        # this part classifies the gene in according to the type of domain and creates a new table with
        # name gene, type domain, pfam id, clas of domain
        
        
        self.input_file_nrgenetics = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOMAIN_LIST_NATREVGENETICS_PROPERTY, True)

        # reading of domain information of dataset 2
        nrgenetics_table = FileParser.make_table(self.input_file_nrgenetics, skip=1)

                


        # extraction of domain domain for each gene
        # extraction of gene id
        domain_column_nrgenetics = TableWrapper.get_column(nrgenetics_table, 2)
        genes_nrgenetics = TableWrapper.get_column(nrgenetics_table, 0)
        
        # This part reads the type of domain for each gene and creates a new columns with the type of domain
        # at the end returns a new table with namegene, domain name, pfam id , class of domain
        #
        # One gene/protein can have more than one domain
        # For each gene the following part checks if the classification of protein domains by comparing the domain with dict_class_domains:
        # 
        # 
        count_classical = 0
        count_nonclassic = 0
        count_unclissified = 0
        count_any_class = 0
        count_no_domain = 0
        new_table_nrgenetics = []
        for i, gene in enumerate(genes_nrgenetics):
            row = []
            row.append(gene)
            row.append(domain_column_nrgenetics[i])
            string_domains = ''
            string_pfamid = ''
            # If the protein hasn't any domain, add static 'no-domains' information
            if domain_column_nrgenetics[i] == '.':
                string_domains += 'No-domains  '
                string_pfamid += '.,'
                count_no_domain += 1
            # If the protein contains some domains checks the class and 
            # akes a string containing the class domains separated by a comma
            else:
                # the domains are separated by a comma
                domains = domain_column_nrgenetics[i].split(',')
                for type_domain in domains:
                    # the domain present also the number information
                    # X-domain[n] for this reason in order to take only the domain name 
                    # the domain is split to '['
                    type_domain = type_domain.split('[')[0]
                    if type_domain in dict_pfam:
                        string_pfamid += dict_pfam[type_domain]+','
                    else:
                        string_pfamid += '-,'
                    print type_domain
                    if type_domain in dict_type_domain:
                        class_domain = dict_type_domain[type_domain][0]
                        if class_domain == 'Classical':
                            string_domains += 'Classical,'
                            count_classical += 1
                        elif class_domain == 'Non-classical':
                            string_domains += 'Non-classical,'
                            count_nonclassic += 1
                        elif class_domain == 'Unclassified':
                            string_domains += 'Unclassified,'
                            count_unclissified += 1
                    elif type_domain not in dict_type_domain:
                        string_domains += 'Other-Domain,'
                        count_any_class += 1                        
                    else:
                        print 'unexpected case', type_domain
            row.append(string_pfamid[0:len(string_pfamid)-1])
            row.append(string_domains[0:len(string_domains)-1])
            new_table_nrgenetics.append(row)
            
            
            
        # print of proteins number for each domain class
           
        print count_classical
        print count_nonclassic
        print count_unclissified
        print count_any_class
        print count_no_domain
            
        
        
        self.path_ouput_file_nrgenetics = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOMAIN_FINAL_TABLE_NATREVGENETICS, True)
        FileWriter.write_table(self.path_ouput_file_nrgenetics, new_table_nrgenetics)
            
            
        # reading of sequences file in order to take the headers 
        # reading of gene and protein id correspondences for different genes
        
        self.file_sequences = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOMAIN_FILE_SEQ_PROPERTY, True)
        self.header = InfoFasta.get_header(self.file_sequences, text=False)
        gene_seq = [ item.split('>')[1].split('|')[0] for item in self.header]
        prot_seq = [ item.split('>')[1].split('|')[2] for item in self.header]
        
        
        # Construction of table containing the RBP protein and the corresponding class domain and pfam id
        count_j = 0
        count_n = 0
        final_table = []
        title = ['gene id', 'type_domain', 'pfam id', 'class domain', 'prot id']
        final_table.append(title)
        # table construction 
        for n, gene in enumerate(gene_seq):
            print n+1
            print gene
            # if gene is in jprot and in nrgenetics
            if gene in genes_jprot and gene in genes_nrgenetics:
                count_n += 1
                ind_gene = genes_nrgenetics.index(gene)
                row = new_table_nrgenetics[ind_gene]
                row.append(prot_seq[n])
            # if gene is in jprot and not in nrgenetics
            elif gene in genes_jprot and gene not in genes_nrgenetics:
                count_j += 1
                ind_gene = genes_jprot.index(gene)
                row = new_table_jprot[ind_gene]
                row.append(prot_seq[n])
            # if gene not in jprot and gene in nrgenetics
            elif gene not in genes_jprot and gene in genes_nrgenetics:
                count_n += 1
                ind_gene = genes_nrgenetics.index(gene)
                row = new_table_nrgenetics[ind_gene]
                row.append(prot_seq[n])
            # if the gene is not in both dataset: Error
            else:
                print 'Error',  gene, prot_seq[n]
            final_table.append(row)
            
        
        
        sort_table = TableWrapper.inv_column(final_table, 0, 4)
        
        sort_table = TableWrapper.inv_column(sort_table, 1,4)
        sort_table = TableWrapper.inv_column(sort_table, 2,4)
        sort_table = TableWrapper.inv_column(sort_table, 3,4)
        
        
        
        self.final_file_table = self.path_home + PropertyManager.get_instance().get_property( DataConstants.DOMAIN_FINALE_TABLE_RBP_DATASET_PROPERTY, True)
        FileWriter.write_table(self.file_final_table, sort_table)
     

    
        # Counting the number of protein with Classical domain, Non-classical, unclassified or combinations of this class
        prot_name = TableWrapper.get_column(sort_table, 1, 1)
        class_domain = TableWrapper.get_column(sort_table, 4, 1)
        classical = []
        nonclassical = []
        unclissified = []
        anyclass = []
        nodomain = []

        # there are several combinations of class domain in a protein
        # the list creation in according to class domain has been performed by the filters of spreadsheet (see documentation)
        for n, domain in enumerate(class_domain):
            diff_domain = domain.split(',')
            for item in diff_domain:
                if item == 'Classical':
                    classical.append(prot_name[n])
                elif item == 'Non-classical':
                    nonclassical.append(prot_name[n])
                elif item == 'Unclassified':
                    unclissified.append(prot_name[n])
                elif item == 'Other-Domain':
                    anyclass.append(prot_name[n])
                elif item == 'No-domains':
                    nodomain.append(prot_name[n])
                else:
                    # many others
                    pass
                
                    
                    

