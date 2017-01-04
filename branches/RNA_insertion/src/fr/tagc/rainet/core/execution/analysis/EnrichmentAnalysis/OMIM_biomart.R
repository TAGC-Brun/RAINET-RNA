#!/usr/bin/R
# 04-Jan-2017 Diogo Ribeiro
# Script to access Ensembl BioMart, retrieve and write OMIM-related fields

# Download / load required libraries
#########################################################################

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

library("biomaRt")

# Setting up constant variables
##########################################################################

#Note: as of 19-Jan-2016, it is not possible to access BioMart through the www.biomart.org/ , only access is through www.ensembl.org
ENSEMBL_URL = "ensembl.org" #this is link to current Ensembl version 87
ENSEMBL_BIOMART = "ENSEMBL_MART_ENSEMBL"
BIOMART_DATASET = "hsapiens_gene_ensembl"

BIOMART_CENTRAL_ATTRIBUTE = "ensembl_transcript_id"

WANTED_ATTRIBUTES =  c(
#  "ensembl_peptide_id",
  "mim_morbid",
#    "hgnc_symbol",
  "uniprot_swissprot"
#  "mim_gene_accession",
  )

OUTPUT_ATTRIBUTES = "/home/diogo/Documents/RAINET_data/OMIM/OMIM_biomart.txt"

BATCH_SIZE = 500 #500 is the advised maximum limit of filter items from BioMart

# Establishing server connection
##########################################################################

ensembl = useMart(host= ENSEMBL_URL, biomart= ENSEMBL_BIOMART, dataset= BIOMART_DATASET)
paste("Accessing host:",ensembl@host)

# Getting all transcript IDs
##########################################################################

allAtt = getBM(attributes=WANTED_ATTRIBUTES,mart = ensembl)
paste("Total number of entries:",nrow(allAtt))

# write to output
##########################################################################
header = capture.output(cat(WANTED_ATTRIBUTES,sep ="\t") )
write.table(header,file = OUTPUT_ATTRIBUTES,quote = FALSE,sep="\t", append = TRUE,  row.names = FALSE, col.names = FALSE)
write.table(allAtt,file = OUTPUT_ATTRIBUTES, quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE, append = TRUE)

print ("FINISHED!")

