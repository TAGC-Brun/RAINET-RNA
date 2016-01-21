# =========================================================================
# This class allows to access to Ensembl (release 75) through URL connection
# =========================================================================

from core.util.web.Web import Web
from core.util.parsing.FileParser import FileParser
from core.util.file.FileUtils import FileUtils
from core.util.log.Logger import Logger
from core.util import Constants


class Ensembl():
    
    def __init__(self):
        
        self.server = Constants.SERVER_RELEASE_75
        self.end = Constants.END_URL_ENSEMBL
        self.dict_query = Constants.TYPE_QUERY_ENSEMBL
        self.section = Constants.SECTION_ENSEMBL
    
    # connection
    # --------------
    #    
    # This method connects to Ensembl
    #
    # Arguments:
    #    - query:  can be a gene or trancript identifier
    # Returns:
    #    - the text page obtained from the request to Ensembl
    #
    def connection(self, query):
        
        query_url = self.dict_query[query[0:4]] + query      
        url = self.server + self.section + query_url + self.end
        # Call the Web.openurl method that can be access to an URL
        text_page = Web.openurl(url)
        
        return text_page
    
    
    # get_sequence
    # ----------------
    # 
    # This method gets all isoform sequences of a gene when:
    # query = gene 
    # query = transcript and num = 'all'
    # 
    # This method gets just one sequence when:
    # query = transcript and num = 'one'
    #
    # Returns:
    #    - a string in fasta format
    
    @staticmethod 
    def get_sequence(query, num='all'): 
        
        E = Ensembl()
        text_page = E.connection(query)
        #
        # The request may return three situation:
        # 1 The gene is not anymore available
        # 2 The gene is a pseudogene therefore any sequence correspond to it
        # 3 The gene have at least one sequence
        if 'Sorry, you have entered an invalid URL' in text_page:
            Logger.get_instance().info(" Ensembl : This genes is not anymore available "  + str(query))
            return query + ' No available'
        elif '>' not in text_page:
            Logger.get_instance().info(" Ensembl : This gene is a pseudogene " + str(query))
            return query + ' pseudogene'
        else:
            # The \r chars are deleted
            text_page = text_page.replace('\r', '') 
            # 
            # Split text page in order to have each isoform as item of list
            #
            list_isoform =text_page.split('>')[1:]
            all_isoform = []
            for item in list_isoform:
                isoform = ''
                item = item.strip()
                row = item.split('\n')
                header = row[0]
                # Definition of identifier in order to obtain a final header like
                # >ENSG|ENST|ENSP
                gene =  header[header.index('ENSG'): header.index('ENSG') +15]
                transcript = header[header.index('ENST'): header.index('ENST') +15]
                protein = header[header.index('ENSP'): header.index('ENSP') +15]
                new_header = '>' + gene + '|' + transcript + '|' + protein
                # Replace the old header with new one
                # The rows (of one isoform) are linked to form a string
                row[0]= new_header
                isoform = '\n'.join(row)
                all_isoform.append(isoform)
            # 
            # All isoform (of gene) are linked to form a string
            #
            string = '\n'.join(all_isoform)
        
            # Different return cases
            #
            if query[0:4] == 'ENSG' and num == 'all':
                return string
            elif query[0:4] == 'ENST' and num =='one':            
                seq = [item for item in all_isoform if query in item]
                return seq[0]
            elif query[0:4] == 'ENST' and num == 'all':
                return string
            elif query[0:4] == 'ENSG' and num == 'one':
                Logger.get_instance().warning(" This combination of args doesn't make sense\n try again!")
            else:
                Logger.get_instance().warning(" The args are wrong, try again!")
                
    
    # get_header
    # ----------------------
    #    
    # This method get the header of gene or transcript 
    # Arguments:
    #    -query: gene or transcript
    #    -amount: 
    #         True --> return all header of gene 
    #         False--> return only the header of the transcript gived as input
    #
    @staticmethod
    def get_header(query, amount=True):
        
        if amount==True:
            E = Ensembl()
            fasta = E.get_sequence(query)
            lines = fasta.split('\n')
            all_header = [row for row in lines if row[0:1]=='>']
            headers = '\n'.join(all_header)
            return headers
        else:
            E = Ensembl()
            fasta = E.get_sequence(query, num='one')
            header = fasta[0:48]
            return header
        
    # list_get_seq
    # ------------------
    # 
    # This method get the sequence of a gene or transcript list
    # Arguments:
    #    - path_input: file containing gene or transcript list
    #    - type_query: 
    #         1 -->all or 2-->one in order to get the sequence
    # optional Arguments
    #    - path_protein_list: protein list of which you want the sequences
    #    - path_output: output file name
    #
    # Return
    #    - file containing the sequences in fasta format
    #
    @staticmethod
    def list_get_seq(path_input, type_query, path_protein_list=None, path_output=None):
        
        # the input file is read
        list_item = FileParser.read_file(path_input)
        dict_query = {1: 'all', 2: 'one'}
        count_duplicate_genes = 0 
        all_seqs =''
        prot_seq = []
        # For each gene in list the sequences are downloaded 
        for i, item in enumerate(list_item):
            Logger.get_instance().info(str(i+1) + ' Extraction of gene sequence(s) : ' + item)
            fasta_seq = Ensembl.get_sequence(item, dict_query[type_query])
            if fasta_seq == item + ' No available':
                pass
            elif fasta_seq == item + ' pseudogene':
                pass
            
            # If the gene have a sequence the output is memorized in seqs
            else:
                seqs = fasta_seq
                seqs = seqs + '\n'
                if path_protein_list == None:
                    pass
                
                # if path_protein_list is different to None
                # Among the isoform of gene will be get only that is present in list_protein
                else:
                    list_protein = FileParser.read_file(path_protein_list)
                    prot_genes = seqs.split('>')
                    protein_seq = ['>'+fasta for fasta in prot_genes if fasta[32:47] in list_protein]
                #    
                # if path_output == None the information are stored in list o string
                #
                if path_output == None:
                    if path_protein_list == None:
                        all_seqs+=seqs
                    else:
                        prot_seq.append(protein_seq)
                #
                # if path_output != None
                # the information will be appended in a file   
                else:
                    file_out = FileUtils.open_text_a(path_output)
                    if path_protein_list == None:
                        file_out.write(seqs)                    
                    else:
                        if protein_seq == []:
                            count_duplicate_genes += 1
                            Logger.get_instance().info(" Number of duplicated genes: " + str(count_duplicate_genes))
                            Logger.get_instance().info(" The gene duplicated is: " + str(item) + '|' + str(list_protein[i]))
                        else:
                            file_out.write(protein_seq[0])
                    
        
        # return information like string or list
                       
        if path_output == None:
            if path_protein_list == None:
                all_seqs+=seqs
                return seqs
            else:
                return protein_seq
        else:       
            file_out.close()
            

        
            
        
        
    
        
        
        
        
        
        
        
        
