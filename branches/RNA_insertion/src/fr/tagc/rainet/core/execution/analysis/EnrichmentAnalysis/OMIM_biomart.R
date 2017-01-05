#!/usr/bin/R
# 04-Jan-2017 Diogo Ribeiro
# Script to access Ensembl BioMart, retrieve and write OMIM-related fields

# Download / load required libraries
#########################################################################

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

library("biomaRt")

# Get input arguments
#########################################################################

# Get the arguments from the launch command line
args <- commandArgs(TRUE)

# Test if we have enough arguments
if( length(args) != 1){
  stop("Rscript: Bad argument number")
}

outputFile = args[1]  #"/home/diogo/Documents/RAINET_data/OMIM/OMIM_biomart.txt"

# Setting up constant variables
##########################################################################

ENSEMBL_URL = "ensembl.org" #this is link to current Ensembl version 87
ENSEMBL_BIOMART = "ENSEMBL_MART_ENSEMBL"
BIOMART_DATASET = "hsapiens_gene_ensembl"

WANTED_ATTRIBUTES =  c(
  "ensembl_gene_id",
  "mim_morbid",
  "hgnc_symbol",
  "uniprot_swissprot"
#  "mim_gene_accession",
  )

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
write.table(header,file = outputFile,quote = FALSE,sep="\t", append = TRUE,  row.names = FALSE, col.names = FALSE)
write.table(allAtt,file = outputFile, quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE, append = TRUE)

print ("FINISHED!")

