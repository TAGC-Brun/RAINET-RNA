#!/usr/bin/env python3

#===============================================================================
# 01-Fev-2016 Diogo Ribeiro
# Script to create test files for Rainet project
# Objective is to be able to run Rainet code under controlled and tested conditions
#===============================================================================

import os
import re
import subprocess
import warnings

from fr.tagc.rainet.core.data import DataConstants
from unittest.main import TestProgram
from fr.tagc.rainet.core.util.exception.RainetException import RainetException

#===============================================================================
PATHS = [] # will store all the file paths for creating a testing .ini file
OUTPUT_FOLDER = "/home/diogo/workspace/tagc-rainet-RNA/test/unittest/testInput/" #note that output of this script will be used as input for testing
INPUT_FOLDER = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/"
SAMPLE_NUMBER = 100
#===============================================================================


# [PROTEINS]
# PROTEIN_UNIPROT_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/human_uniprot_protein_list.txt
# PROTEIN_CROSSREFERENCES = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/HUMAN_9606_idmapping.dat
# PROTEIN_ISOFORMS = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/UP000005640_9606_additional.fasta
# PROTEIN_DOMAIN_SMART = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/all_protein_domains_smart.txt
# PROTEIN_DOMAIN_PFAM = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/Pfam-A.clans.tsv
# 
# [GENE_ONTONLOGY]
# GENE_ONTOLOGY_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/go-basic.obo
# GENE_ONTOLOGY_ANNOTATION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/gene_association.goa_human
# 
# [KEGG_PATHWAY]
# KEGG_PATHWAY_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/human_kegg_pathway_definitions.txt
# KEGG_PATHWAY_ANNOTATION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/human_kegg_pathway_annotations.txt
# 
# [REACTOME PATHWAY]
# REACTOME_PATHWAY_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/all_reactome_pathway_definitions.txt
# REACTOME_PATHWAY_ANNOTATION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/all_reactome_pathway_annotations.txt
# 
# [INTERACTOME]
# INTERACTOME_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/human.pairmap
# INTERACTOME_NETWORK_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/human.binary.nr0.95.connected.gr
# INTERACTOME_NETWORK_PARTITION_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/human.binary.nr0.95.connected.clas
# INTERACTOME_NETWORK_PARTITION_ANNOTATION =  /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/human.binary.nr0.95.connected.fm
# INTERACTOME_NETWORK_REDUNDANCY_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/human_0.95.blastmap
# 
# [RNA]
# RNA_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/RNA/RNA_ATTRIBUTES.tsv
# RNA_CROSS_REFERENCE = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/RNA/RNA_XREF_ATTRIBUTES.tsv
# 
# [PROTEIN_RNA_INTERACTION]
# PROTEIN_RNA_INTERACTION_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/PRI/outFile_T-T.2825-7004.sorted.txt


#===============================================================================
# Approach, pick 100 proteins and 100 RNAs randomly from initial input files, grep data for them in all other files
#===============================================================================

#===============================================================================
# Files that do not need to be newly created:
# 
# all_protein_domains_smart.txt
# Pfam-A.clans.tsv
# go-basic.obo
# human_kegg_pathway_definitions.txt
# human_kegg_pathway_annotations.txt
# all_reactome_pathway_definitions.txt
# all_reactome_pathway_annotations.txt
#===============================================================================

def runProcess(cmd):
    """General function to run external command and return the return code, and standard outputs"""

    run = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    run.wait()
  
    stderrText = run.stderr.read().decode("UTF-8").strip()
    
    if len(stderrText) > 0:
        warnings.warn("STDERR:\n"+stderrText)

    if run.returncode != 0:
        warnings.warn("ERROR: Return code:"+str(run.returncode)+"\t"+cmd)
  
    stdoutText = run.stdout.read().decode("UTF-8").strip()

    return stdoutText


def makeProteinUniprotDefinition(sample_number):

    #===============================================================================
    # File: human_uniprot_protein_list.txt ; PROTEIN_UNIPROT_DEFINITION
    #===============================================================================
    
    #Example lines
    #Entry   Entry name      Protein names   Gene names  (primary )  Gene names  (synonym )  Organism        Length  Fragment        Cross-reference (Pfam)  Cross-reference (SMART)
    #P31946  1433B_HUMAN     14-3-3 protein beta/alpha (Protein 1054) (Protein kinase C inhibitor protein 1) (KCIP-1) [Cleaved into: 14-3-3 protein beta/alpha, N-terminally processed]      YWHAB           Homo sapiens (Human)    246             PF00244;        SM00101;
    
    #Columns
    #PROTEIN_HEADERS = ["Entry", "Entry name", "Protein names", "Gene names  (primary )", "Gene names  (synonym )", "Organism", "Length", "Fragment", "Cross-reference (Pfam)", "Cross-reference (SMART)"]
    #PROTEIN_PARAMS = ["Entry", "Entry name", "Protein names", "Gene names  (primary )", "Gene names  (synonym )", "Organism", "Length", "Fragment", "Cross-reference (Pfam)", "Cross-reference (SMART)"]

    print ("makeProteinUniprotDefinition..")

    testingProteins = set() # will store the initial sample of proteins

    inFile = INPUT_FOLDER+"human_uniprot_protein_list.txt"
    outFile = OUTPUT_FOLDER+"human_uniprot_protein_list.txt" #"protein_uniprot.txt"

    outHandler = open(outFile,"w")
    outHandler.write("\t".join(DataConstants.PROTEIN_HEADERS)+"\n")
    outHandler.close()

    cmd = "shuf -n %s %s >> %s"  % (sample_number,inFile,outFile)
    runProcess(cmd)



def makeProteinCrossReferences(testingProteins):

    #===============================================================================
    # File: HUMAN_9606_idmapping.dat ; PROTEIN_CROSSREFERENCES
    #===============================================================================
    
    # Example
    # P31946  UniProtKB-ID    1433B_HUMAN
    # P31946  GI      4507949
    
    proteinXref = {} #key -> xref, val -> main protein id
    
    print("makeProteinCrossReferences..")
    
    inFile = INPUT_FOLDER+"HUMAN_9606_idmapping.dat"
    outFile = OUTPUT_FOLDER+"HUMAN_9606_idmapping.dat" #"protein_crossreferences.txt"
     
    outHandler = open(outFile,"w")
     
    for protein in testingProteins:
        cmd = "awk '$1==%s' %s" % ('"'+protein+'"',inFile)
#        cmd = "grep %s %s" % ('"'+protein+'"',inFile) #this is faster but may bring other entries that do not correspond to our chosen proteins
        out = runProcess(cmd)
        ids = out.strip().split("\n")
        for items in ids:
            id = items.split("\t")[2]
            proteinXref[id] = protein
        outHandler.write(out+"\n")
     
    outHandler.close()

    return proteinXref


def makeProteinIsoforms(testingProteins):
    
    #===============================================================================
    # File: UP000005640_9606_additional.fasta ; PROTEIN_ISOFORMS
    #===============================================================================

    # Example
    # >tr|A0A024QZ18|A0A024QZ18_HUMAN Isoform of Q9H6Y5, PDZ domain-containing protein MAGIX OS=Homo sapiens GN=MAGIX PE=1 SV=1
    # MPLLWITGPRYHLILLSEASCLRANYVHLCPLFQHRWLETCNAPPQLIQGKARSAPKPSQ

    #Dependency: fastaq, from https://github.com/sanger-pathogens/Fastaq

    print("makeProteinIsoforms..")

    inFile = INPUT_FOLDER+"UP000005640_9606_additional.fasta"
    tmpFile1 = OUTPUT_FOLDER+"makeProteinIsoforms.tmp1"
    tmpFile2 = OUTPUT_FOLDER+"makeProteinIsoforms.tmp2"
    outFile = OUTPUT_FOLDER+"UP000005640_9606_additional.fasta" #"protein_isoforms.fasta"
    
    cmd = "fastaq get_ids %s %s" % (inFile,tmpFile1)
    runProcess(cmd)

    tmpHandler = open(tmpFile2,"w")
    
    for protein in testingProteins:
        cmd = "grep %s %s" % (protein,tmpFile1) # grep returns code 1 if not finding anything..
        tmpHandler.write(runProcess(cmd))

    tmpHandler.close()

    cmd = "fastaq filter --ids_file %s %s %s" % (tmpFile2,inFile,outFile)
    runProcess(cmd)

    os.remove(tmpFile1)
    os.remove(tmpFile2)
   

def makeGeneOntologyAnnotation(testingProteins):
    
    #===============================================================================
    # File: gene_association.goa_human ; GENE_ONTOLOGY_ANNOTATION
    #===============================================================================

    # Example
    #UniProtKB       A0A024QZ42      PDCD6           GO:0005509      GO_REF:0000002  IEA     InterPro:IPR002048      F       HCG1985580, isoform CRA_c       A0A024QZ42_HUMAN|PDCD6|hCG_1985580      protein taxon:9606      20160102        InterPro      

    print("makeGeneOntologyAnnotation..")

    inFile = INPUT_FOLDER+"gene_association.goa_human"
    outFile = OUTPUT_FOLDER+"gene_association.goa_human" #"gene_ontology_annotation.goa_human"

    inHandler = open(inFile,"r")    
    outHandler = open(outFile,"w")

    #writing header    
    for line in inHandler:
        if line.startswith(DataConstants.PROTEIN_GO_ANNOTATION_COMMENT_CHAR):
            outHandler.write(line)
        else:
            spl = line.strip().split("\t")
            if spl[1] in testingProteins:
                outHandler.write(line)
     
#     for protein in testingProteins:
#         cmd = "awk '$2==%s' %s" % ('"'+protein+'"',inFile)
#         outHandler.write(runProcess(cmd))
     
    outHandler.close()


def makeInteractomeDefinition(testingProteins,testingXref):
    
    #===============================================================================
    # File: human.pairmap ; INTERACTOME_DEFINITION
    #===============================================================================

    # Example
    #     entrez gene/locuslink:5764      PTN_HUMAN       entrez gene/locuslink:6745      SSRA_HUMAN      MI:0018(two hybrid)     PUBMED:16169070 BIOGRID  "MI:0407"(direct interaction)  score:3.0
    #     entrez gene/locuslink:10865     ARI5A_HUMAN     entrez gene/locuslink:2100      ESR2_HUMAN      MI:0018(two hybrid)     PUBMED:15941852 BIOGRID  "MI:0407"(direct interaction)  -

    print ("makeInteractomeDefinition..")

    inFile = INPUT_FOLDER+"human.pairmap"
    outFile = OUTPUT_FOLDER+"human.pairmap" #"interactome_definition.pairmap"

    inHandler = open(inFile,"r")    
    outHandler = open(outFile,"w")

    interactionsDict = {}
    for line in inHandler:
        spl = line.strip().split("\t")
        tag = spl[0]+"|"+spl[2]
        if tag not in interactionsDict:
            interactionsDict[tag] = []
        interactionsDict[tag].append(line)

    #because interaction file seems to use all kinds of IDs, I should be using xrefs for getting the maximum protein-protein interactions
    testingXref = set(testingProteins).union(set(proteinXref.keys()) )

    for protein1 in testingXref:
        for protein2 in testingXref:
            pair1 = protein1+"|"+protein2
            pair2 = protein2+"|"+protein1
            if pair1 in interactionsDict:
                for p in interactionsDict[pair1]:
                    outHandler.write(p)
            if pair2 in interactionsDict:
                for p in interactionsDict[pair2]:
                    outHandler.write(p)
            
    outHandler.close()


def makeInteractomeNetworkDefinition(testingProteins,proteinXref):
    
    #===============================================================================
    # File: human.binary.nr0.95.connected.gr ; INTERACTOME_NETWORK_DEFINITION
    #===============================================================================

    # Example
    #     PRRT3_HUMAN     TMM17_HUMAN
    #     HES6_HUMAN      TLE1_HUMAN

    print ("makeInteractomeNetworkDefinition..")

    inFile = INPUT_FOLDER+"human.binary.nr0.95.connected.gr"
    outFile = OUTPUT_FOLDER+"human.binary.nr0.95.connected.gr" #"interactome_network_definition.pairmap"

    inHandler = open(inFile,"r")    
    outHandler = open(outFile,"w")

    interactionsDict = {}
    for line in inHandler:
        spl = line.strip().split("\t")
        tag = spl[0]+"|"+spl[1]
        if tag not in interactionsDict:
            interactionsDict[tag] = []
        interactionsDict[tag].append(line)

    #because interaction file seems to use all kinds of IDs, I should be using xrefs for getting the maximum protein-protein interactions
    testingXref = set(testingProteins).union(set(proteinXref.keys()) )

    for protein1 in testingXref:
        for protein2 in testingXref:
            pair1 = protein1+"|"+protein2
            pair2 = protein2+"|"+protein1
            if pair1 in interactionsDict:
                for p in interactionsDict[pair1]:
                    outHandler.write(p)
            if pair2 in interactionsDict:
                for p in interactionsDict[pair2]:
                    outHandler.write(p)
            
    outHandler.close()



#===============================================================================
# Client
#===============================================================================

makeProteinUniprotDefinition(SAMPLE_NUMBER)

#getter
testingProteins = []
proteinListHandler = open(OUTPUT_FOLDER+"protein_list/testing_proteins_list.txt","w")
with open(OUTPUT_FOLDER+"human_uniprot_protein_list.txt","r") as f:
    f.readline() #skip header
    for line in f:
        prot = line.strip().split("\t")[0]
        testingProteins.append(prot)
        proteinListHandler.write(prot+"\n")
proteinListHandler.close()
assert(len(testingProteins) == SAMPLE_NUMBER),"Number of sampled proteins should be equal to given SAMPLE_NUMBER"

proteinXref = makeProteinCrossReferences(testingProteins)
makeProteinIsoforms(testingProteins)
makeGeneOntologyAnnotation(testingProteins)
makeInteractomeDefinition(testingProteins,proteinXref)
makeInteractomeNetworkDefinition(testingProteins,proteinXref)

print ("FINISHED!")

# #02-Fev-2016
# #"makeTestFiles.py"
# #I tried to create a smaller version of the database, by making new files, with information only on a small subset of proteins.
# #This turned out to be quite a major undertaking. I've only make the following files:
# human_uniprot_protein_list.txt
# HUMAN_9606_idmapping.dat
# UP000005640_9606_additional.fasta
# gene_association.goa_human
# human.pairmap
# human.binary.nr0.95.connected.gr
# #all other input files were copied straight from Lionel
# 
# #Insertion of all this data took  9 minutes 51 seconds







