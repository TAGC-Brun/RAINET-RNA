#!/usr/bin/env python3

#===============================================================================
# 01-Fev-2016 Diogo Ribeiro
# Script to create test files for Rainet project
#
# Objective is to be able to run Rainet code under controlled and tested conditions
# Approach: pick X proteins and X RNAs randomly from initial input files (Protein models / RNA models), 
# grep/awk data associated to them in all other files
#===============================================================================

import os
import sys
import re
import subprocess
import warnings

from fr.tagc.rainet.core.data import DataConstants

if sys.version_info.major < 3:
   warnings.warn("This script uses only Python 3, check your version. Exiting..")
   sys.exit(1)

#===============================================================================
# Constants
#
#===============================================================================

# Store all the files that will be copied as they are originally (files that we do not make smaller for testing database)
REMAINING_FILES = ["all_protein_domains_smart.txt",
                   "Pfam-A.clans.tsv","go-basic.obo",
                   "human_kegg_pathway_definitions.txt",
                   "human_kegg_pathway_annotations.txt",
                   "all_reactome_pathway_definitions.txt",
                   "human.binary.nr0.95.connected.noself.gr",
                   "human.binary.nr0.95.connected.noself.clas",
                   "human.binary.nr0.95.connected.noself.fm",
                   "human_0.95.blastmap"]

REMAINING_FILES_RNA = ["list_interating_RNAs.txt","list_interating_proteins.txt"]

#note that output of this script will be used as input for testing
OUTPUT_FOLDER = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/" 

INPUT_FOLDER = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/input_data/"

#number of proteins and RNAs for testing
SAMPLE_NUMBER = 200

#Use 1 if wanting to recreate random list of proteins/RNAs
FORCE_CREATION = 0

#===============================================================================

if not os.path.exists(OUTPUT_FOLDER+"/PROTEIN"): os.mkdir(OUTPUT_FOLDER+"/PROTEIN")
if not os.path.exists(OUTPUT_FOLDER+"/RNA"): os.mkdir(OUTPUT_FOLDER+"/RNA")

#===============================================================================
# Copy of file used for full Rainet DB insertion. We want to make smaller file versions for some of these.
#
#===============================================================================
# [PROTEINS]
# PROTEIN_UNIPROT_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/human_uniprot_protein_list.txt
# PROTEIN_CROSSREFERENCES = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/HUMAN_9606_idmapping.dat
# PROTEIN_ISOFORMS = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/UP000005640_9606_additional.fasta
# PROTEIN_DOMAIN_SMART = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/all_protein_domains_smart.txt
# PROTEIN_DOMAIN_PFAM = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/Pfam-A.clans.tsv
# 
# [GENE_ONTONLOGY]
# GENE_ONTOLOGY_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/go-basic.obo
# GENE_ONTOLOGY_ANNOTATION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/gene_association.goa_human
# 
# [KEGG_PATHWAY]
# KEGG_PATHWAY_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/human_kegg_pathway_definitions.txt
# KEGG_PATHWAY_ANNOTATION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/human_kegg_pathway_annotations.txt
# 
# [REACTOME PATHWAY]
# REACTOME_PATHWAY_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/all_reactome_pathway_definitions.txt
# REACTOME_PATHWAY_ANNOTATION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/all_reactome_pathway_annotations.txt
# 
# [INTERACTOME]
# INTERACTOME_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/human.pairmap
# INTERACTOME_NETWORK_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/human.binary.nr0.95.connected.noself.gr
# INTERACTOME_NETWORK_PARTITION_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/human.binary.nr0.95.connected.noself.clas
# INTERACTOME_NETWORK_PARTITION_ANNOTATION =  /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/human.binary.nr0.95.connected.noself.fm
# INTERACTOME_NETWORK_REDUNDANCY_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/PROTEIN/human_0.95.blastmap
# 
# [RNA]
# RNA_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/RNA/RNA_ATTRIBUTES_ensembl82.tsv
# RNA_CROSS_REFERENCE = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/RNA/RNA_XREF_ATTRIBUTES_ensembl82.tsv
# 
# [PROTEIN_RNA_INTERACTION]
# PROTEIN_RNA_INTERACTION_DEFINITION = /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/RNA/outFile_T-T.2825-7004.sorted.txt

#===============================================================================
# Files that do not need to be newly created:
#===============================================================================
#
# all_protein_domains_smart.txt
# Pfam-A.clans.tsv
# go-basic.obo
# human_kegg_pathway_definitions.txt
# human_kegg_pathway_annotations.txt
# all_reactome_pathway_definitions.txt
# all_reactome_pathway_annotations.txt
# human.binary.gr
# human.binary.nr0.95.connected.clas
# human.binary.nr0.95.connected.fm
# human_0.95.blastmap
#===============================================================================


# #
# General function to run command using subprocess, return only standard output, but take care of standard error or run failures
def run_process(cmd):
    """General function to run external command and return the return code, and standard outputs"""

    run = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    run.wait()
  
    stderrText = run.stderr.read().decode("UTF-8").strip()
    
    if len(stderrText) > 0:
        warnings.warn("STDERR:\n"+stderrText)

    if run.returncode != 0:
        warnings.warn("ERROR: Return code:"+str(run.returncode)+"\t"+cmd)
  
    stdoutText = run.stdout.read().decode("UTF-8").strip()

#     if stdoutText == None or len(stdoutText) == 0:
#         warnings.warn("Output line is empty :"+str(cmd))

    return stdoutText


# #
# Getter for making file with list of proteins and assert correctness
# By doing this I can rerun without changing the initial set of chosen proteins
def list_getter(sample_number, in_file):

    listOfItems = set()
    
    with open(in_file,"r") as f:
        for line in f:
            prot = line.strip()
            listOfItems.add(prot) #.append(prot)

    assert(len(listOfItems) == SAMPLE_NUMBER),"Number of items sampled should be equal to given SAMPLE_NUMBER"

    return listOfItems


def copy_files(list_of_files,in_folder, out_folder):

    for file in list_of_files:
        cmd = "cp %s %s" % (in_folder+file,out_folder)
        run_process(cmd)


def make_sample(sample_number, inFile, outFolder, outFile):

    #First, run shuffling command to create initial file
    print ("make_sample_file..")

    if not os.path.exists(outFolder):
        os.mkdir(outFolder)

    cmd = "shuf -n %s %s | cut -f1 > %s"  % (sample_number,inFile,outFolder+"/"+outFile)
    run_process(cmd)


def make_protein_uniprot_definition(testingProteins):

    #===============================================================================
    # File: human_uniprot_protein_list.txt ; PROTEIN_UNIPROT_DEFINITION
    #===============================================================================
    
    #Example lines
    #Entry   Entry name      Protein names   Gene names  (primary )  Gene names  (synonym )  Organism        Length  Fragment        Cross-reference (Pfam)  Cross-reference (SMART)
    #P31946  1433B_HUMAN     14-3-3 protein beta/alpha (Protein 1054) (Protein kinase C inhibitor protein 1) (KCIP-1) [Cleaved into: 14-3-3 protein beta/alpha, N-terminally processed]      YWHAB           Homo sapiens (Human)    246             PF00244;        SM00101;
    
    #Columns
    #PROTEIN_HEADERS = ["Entry", "Entry name", "Protein names", "Gene names  (primary )", "Gene names  (synonym )", "Organism", "Length", "Fragment", "Cross-reference (Pfam)", "Cross-reference (SMART)"]
    #PROTEIN_PARAMS = ["Entry", "Entry name", "Protein names", "Gene names  (primary )", "Gene names  (synonym )", "Organism", "Length", "Fragment", "Cross-reference (Pfam)", "Cross-reference (SMART)"]


    print ("make_protein_uniprot_definition..")

    inFile = INPUT_FOLDER+"PROTEIN/human_uniprot_protein_list.txt"
    outFile = OUTPUT_FOLDER+"PROTEIN/human_uniprot_protein_list.txt" #"protein_uniprot.txt"

    inHandler = open(inFile,"r")
    outHandler = open(outFile,"w")

    outHandler.write(inHandler.readline()) #write header

    for line in inHandler:
        prot = line.split("\t")[0]
        if prot in testingProteins:
            outHandler.write(line)
        
    outHandler.close()


def make_protein_cross_references(testing_proteins):

    #===============================================================================
    # File: HUMAN_9606_idmapping.dat ; PROTEIN_CROSSREFERENCES
    #===============================================================================
    
    # Example
    # P31946  UniProtKB-ID    1433B_HUMAN
    # P31946  GI      4507949
    
    proteinXref = {} #key -> xref, val -> main protein id
    
    print("make_protein_cross_references..")
    
    inFile = INPUT_FOLDER+"PROTEIN/HUMAN_9606_idmapping.dat"
    outFile = OUTPUT_FOLDER+"PROTEIN/HUMAN_9606_idmapping.dat" #"protein_crossreferences.txt"
     
    outHandler = open(outFile,"w")
     
    for protein in testing_proteins:
        cmd = "awk '$1==%s' %s" % ('"'+protein+'"',inFile)
#        cmd = "grep %s %s" % ('"'+protein+'"',inFile) #this is faster but may bring other entries that do not correspond to our chosen proteins
        out = run_process(cmd)
        ids = out.strip().split("\n")
        for items in ids:
            id = items.split("\t")[2]
            proteinXref[id] = protein
        outHandler.write(out+"\n")
     
    outHandler.close()

    return proteinXref


def make_protein_isoforms(testing_proteins):
    
    #===============================================================================
    # File: UP000005640_9606_additional.fasta ; PROTEIN_ISOFORMS
    #===============================================================================

    # Example
    # >tr|A0A024QZ18|A0A024QZ18_HUMAN Isoform of Q9H6Y5, PDZ domain-containing protein MAGIX OS=Homo sapiens GN=MAGIX PE=1 SV=1
    # MPLLWITGPRYHLILLSEASCLRANYVHLCPLFQHRWLETCNAPPQLIQGKARSAPKPSQ

    #Dependency: fastaq, from https://github.com/sanger-pathogens/Fastaq

    print("make_protein_isoforms..")

    inFile = INPUT_FOLDER+"PROTEIN/UP000005640_9606_additional.fasta"
    tmpFile1 = OUTPUT_FOLDER+"PROTEIN/make_protein_isoforms.tmp1"
    tmpFile2 = OUTPUT_FOLDER+"PROTEIN/make_protein_isoforms.tmp2"
    outFile = OUTPUT_FOLDER+"PROTEIN/UP000005640_9606_additional.fasta" #"protein_isoforms.fasta"
    
    cmd = "fastaq get_ids %s %s" % (inFile,tmpFile1)
    run_process(cmd)

    tmpHandler = open(tmpFile2,"w")
    
    for protein in testing_proteins:
        cmd = "grep %s %s" % (protein,tmpFile1) # grep returns code 1 if not finding anything..
        tmpHandler.write(run_process(cmd))

    tmpHandler.close()

    cmd = "fastaq filter --ids_file %s %s %s" % (tmpFile2,inFile,outFile)
    run_process(cmd)

    os.remove(tmpFile1)
    os.remove(tmpFile2)
   

def make_gene_ontology_annotation(testing_proteins):
    
    #===============================================================================
    # File: gene_association.goa_human ; GENE_ONTOLOGY_ANNOTATION
    #===============================================================================

    # Example
    #UniProtKB       A0A024QZ42      PDCD6           GO:0005509      GO_REF:0000002  IEA     InterPro:IPR002048      F       HCG1985580, isoform CRA_c       A0A024QZ42_HUMAN|PDCD6|hCG_1985580      protein taxon:9606      20160102        InterPro      

    print("make_gene_ontology_annotation..")

    inFile = INPUT_FOLDER+"PROTEIN/gene_association.goa_human"
    outFile = OUTPUT_FOLDER+"PROTEIN/gene_association.goa_human" #"gene_ontology_annotation.goa_human"

    inHandler = open(inFile,"r")    
    outHandler = open(outFile,"w")

    #writing header    
    for line in inHandler:
        if line.startswith(DataConstants.PROTEIN_GO_ANNOTATION_COMMENT_CHAR):
            outHandler.write(line)
        else:
            spl = line.strip().split("\t")
            if spl[1] in testing_proteins:
                outHandler.write(line)
     
#     for protein in testingProteins:
#         cmd = "awk '$2==%s' %s" % ('"'+protein+'"',inFile)
#         outHandler.write(run_process(cmd))
     
    outHandler.close()


def make_reactome_pathway_annotation(testing_proteins):
    
    #===============================================================================
    # File: all_reactome_pathway_annotations.txt ; REACTOME_PATHWAY_ANNOTATION
    #===============================================================================

    # Example
    # Q6PIL0    R-HSA-2029481    http://www.reactome.org/PathwayBrowser/#/R-HSA-2029481    FCGR activation    TAS    Homo sapiens

    print("make_reactome_pathway_annotation..")

    inFile = INPUT_FOLDER+"PROTEIN/all_reactome_pathway_annotations.txt"
    outFile = OUTPUT_FOLDER+"PROTEIN/all_reactome_pathway_annotations.txt"

    inHandler = open(inFile,"r")    
    outHandler = open(outFile,"w")

    #writing header    
    for line in inHandler:
        spl = line.strip().split("\t")
        if spl[0] in testing_proteins:
            outHandler.write(line)
     
    outHandler.close()


def make_interactome_definition(testing_proteins,protein_xref):
    
    #===============================================================================
    # File: human.pairmap ; INTERACTOME_DEFINITION
    #===============================================================================

    # Example
    #     entrez gene/locuslink:5764      PTN_HUMAN       entrez gene/locuslink:6745      SSRA_HUMAN      MI:0018(two hybrid)     PUBMED:16169070 BIOGRID  "MI:0407"(direct interaction)  score:3.0
    #     entrez gene/locuslink:10865     ARI5A_HUMAN     entrez gene/locuslink:2100      ESR2_HUMAN      MI:0018(two hybrid)     PUBMED:15941852 BIOGRID  "MI:0407"(direct interaction)  -

    print ("make_interactome_definition..")

    inFile = INPUT_FOLDER+"PROTEIN/human.pairmap"
    outFile = OUTPUT_FOLDER+"PROTEIN/human.pairmap" #"interactome_definition.pairmap"

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
    testingXref = set(testing_proteins).union(set(protein_xref.keys()) )

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


def make_interactome_network_definition(testing_proteins,protein_xref):
    
    #===============================================================================
    # File: human.binary.nr0.95.connected.noself.gr ; INTERACTOME_NETWORK_DEFINITION
    #===============================================================================

    # Example
    #     PRRT3_HUMAN     TMM17_HUMAN
    #     HES6_HUMAN      TLE1_HUMAN

    print ("make_interactome_network_definition..")

    inFile = INPUT_FOLDER+"PROTEIN/human.binary.nr0.95.connected.noself.gr"
    outFile = OUTPUT_FOLDER+"PROTEIN/human.binary.nr0.95.connected.noself.gr" #"interactome_network_definition.pairmap"

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
    testingXref = set(testing_proteins).union(set(protein_xref.keys()) )

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


def make_RNA_definition(testing_RNAs):

    #===============================================================================
    # File: RNA_ATTRIBUTES.tsv ; RNA_DEFINITION
    #===============================================================================
    
    #Example lines
    #  ensembl_transcript_id   ensembl_gene_id ensembl_peptide_id      transcript_biotype      transcript_length       transcript_source       transcript_status       transcript_tsl  transcript_gencode_basic    transcript_start        transcript_end  strand  chromosome_name percentage_gc_content
    #  ENST00000000233 ENSG00000004059 ENSP00000000233 protein_coding  1103    ensembl_havana  KNOWN   tsl1 (assigned to previous version 8)   GENCODE basic   127588345       127591705       1       7       56.5
    
    # RNA_HEADERS = ["transcript_ID","parent_gene","peptide_ID","transcript_biotype","transcript_length","transcript_source","transcript_status","transcript_tsl","transcript_gencode_basic","transcript_start","transcript_end","transcript_strand","chromosome_name","percentage_GC_content"]
    # RNA_PARAMS = ["transcript_ID","parent_gene","peptide_ID","transcript_biotype","transcript_length","transcript_source","transcript_status","transcript_tsl","transcript_gencode_basic","transcript_start","transcript_end","transcript_strand","chromosome_name","percentage_GC_content"]

    print ("make_RNA_definition..")

    #First, run shuffling command to get initial sample

    inFile = INPUT_FOLDER+"RNA/RNA_ATTRIBUTES_ensembl82.tsv"
    outFile = OUTPUT_FOLDER+"RNA/RNA_ATTRIBUTES_ensembl82.tsv" 

    inHandler = open(inFile,"r")
    outHandler = open(outFile,"w")

    outHandler.write(inHandler.readline()) #write header

    for line in inHandler:
        RNA = line.split("\t")[0]
        if RNA in testing_RNAs:
            outHandler.write(line)

    outHandler.close()


def make_RNA_cross_references(testing_RNAs):

    #===============================================================================
    # File: RNA_XREF_ATTRIBUTES.tsv ; RNA_CROSS_REFERENCE
    #===============================================================================
    
    # Example
    # ENST00000215906 hgnc_id HGNC:30437
    # ENST00000000412 hpa     HPA040445
    
    print("make_RNA_cross_references..")
    
    inFile = INPUT_FOLDER+"RNA/RNA_XREF_ATTRIBUTES_ensembl82.tsv"
    outFile = OUTPUT_FOLDER+"RNA/RNA_XREF_ATTRIBUTES_ensembl82.tsv"
     
    outHandler = open(outFile,"w")
     
    for RNA in testing_RNAs:
        cmd = "awk '$1==%s' %s" % ('"'+RNA+'"',inFile)
        out = run_process(cmd)
        if len(out) > 1:
            outHandler.write(out+"\n")
     
    outHandler.close()


def make_PRI_catRAPID(testing_RNAs, testing_proteins):

    #===============================================================================
    # File: catRAPID_interactions.txt ; PROTEIN_RNA_INTERACTION_CATRAPID_DEFINITION
    #===============================================================================

    print("make_PRI_catRAPID..")
    
    inFile = INPUT_FOLDER+"/RNA/catRAPID_interactions_test.txt"
    outFile = OUTPUT_FOLDER+"/RNA/catRAPID_interactions_test.txt"

    # Example
    # sp|Q96DC8|ECHD3_HUMAN ENST00000579524   -12.33  0.10    0.00
    # sp|P10645|CMGA_HUMAN ENST00000516610    10.66   0.32    0.00
    # protein and rna separated by " ", other values separated by "\t"
    
    outHandler = open(outFile,"w")
           
    patterns = set()
    for RNA in testing_RNAs:
        for protein in testing_proteins:
            key = protein+"_"+RNA
            patterns.add(key)

    inHandler = open(inFile,"r")
    for line in inHandler:
        spl = line.split(" ")
        pair = spl[0].split("|")[1]+"_"+spl[1].split("\t")[0]
        if pair in patterns:
            outHandler.write(line)

    outHandler.close()


def make_tx_expression(testing_RNAs):

    #===============================================================================
    # File: transcript_expression_metrics_no_outliers.tsv ; RNA_TISSUE_EXPRESSION_PROPERTY 
    #===============================================================================

    # Example
    # transcript_id   tissue_name     rpkm_mean       rpkm_std        rpkm_median     coef_variation  max
    # ENST00000461495 Thyroid 0.000   0.000   0.000   0       0.000
    # ENST00000461495 Testis  0.000   0.000   0.000   0       0.000    

    print("make_tx_expression..")
    
    inFile = INPUT_FOLDER+"/RNA/transcript_expression_metrics_no_outliers.tsv"
    outFile = OUTPUT_FOLDER+"/RNA/transcript_expression_metrics_no_outliers.tsv"
     
    outHandler = open(outFile,"w")
      
    inHandler = open(inFile,"r")
    
    # Get header
    outHandler.write(inHandler.readline())

    # uniqueFound = set()
    
    for line in inHandler:
        tx = line.split("\t")[0]
        if tx in testing_RNAs:
            # uniqueFound.add(tx)
            outHandler.write(line)

    # Note: the testingRNAs may be more than the number of transcripts retrieved from GTEx. This reflects difference in Ensembl versions used.

    outHandler.close()


#===============================================================================
# Client
#
#===============================================================================

#===============================================================================
# Proteins and PPI
#===============================================================================
if FORCE_CREATION:
    make_sample(SAMPLE_NUMBER, INPUT_FOLDER+"PROTEIN/human_uniprot_protein_list.txt", OUTPUT_FOLDER+"/sampled_item_list/","testing_proteins_list.txt")
testingProteins = list_getter(SAMPLE_NUMBER,OUTPUT_FOLDER+"/sampled_item_list/testing_proteins_list.txt")
    
make_protein_uniprot_definition(testingProteins)
proteinXref = make_protein_cross_references(testingProteins)
make_protein_isoforms(testingProteins)
make_gene_ontology_annotation(testingProteins)
make_reactome_pathway_annotation(testingProteins)
make_interactome_definition(testingProteins,proteinXref)
make_interactome_network_definition(testingProteins,proteinXref)
     
#===============================================================================
# RNA and PRI
#===============================================================================
if FORCE_CREATION:
    make_sample(SAMPLE_NUMBER, INPUT_FOLDER+"RNA/RNA_ATTRIBUTES.tsv", OUTPUT_FOLDER+"/sampled_item_list/","testing_RNA_list.txt")
testingRNAs = list_getter(SAMPLE_NUMBER,OUTPUT_FOLDER+"/sampled_item_list/testing_RNA_list.txt")
  
make_RNA_definition(testingRNAs)
make_RNA_cross_references(testingRNAs)
make_PRI_catRAPID(testingRNAs,testingProteins)
make_tx_expression(testingRNAs)
 
# Copy remaining (unchanged) files
copy_files(REMAINING_FILES,INPUT_FOLDER+"PROTEIN/",OUTPUT_FOLDER+"PROTEIN/")
copy_files(REMAINING_FILES_RNA,INPUT_FOLDER+"RNA/",OUTPUT_FOLDER+"RNA/")

#===============================================================================

print ("FINISHED!")

