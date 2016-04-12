#!/usr/bin/R
# 19-Jan-2016 Diogo Ribeiro
# Script to access Ensembl BioMart, retrieve and write wanted fields for RAINET project

# Download / load required libraries
#########################################################################

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

library("biomaRt")

# Setting up constant variables
##########################################################################

#Note: as of 19-Jan-2016, it is not possible to access BioMart through the www.biomart.org/ , only access is through www.ensembl.org
ENSEMBL_URL = "sep2015.archive.ensembl.org" #this is link to Ensembl version 82, released in September 2015
ENSEMBL_BIOMART = "ENSEMBL_MART_ENSEMBL"
BIOMART_DATASET = "hsapiens_gene_ensembl"

BIOMART_CENTRAL_ATTRIBUTE = "ensembl_transcript_id"

RNA_ATTRIBUTES = c(
  "ensembl_transcript_id",
  "ensembl_gene_id",
  "ensembl_peptide_id",
  "transcript_biotype",
  "transcript_length",
  "transcript_source",
  "transcript_status",
  "transcript_tsl",
  "transcript_gencode_basic",
  "transcript_start",
  "transcript_end",
  "strand",
  "chromosome_name",
  "percentage_gc_content",
  "description"
)

RNA_CROSS_REFERENCE_ATTRIBUTES =  c(
  "hgnc_id",
  "hpa",
  "ucsc",
  "embl",
  "ccds",
  "uniprot_swissprot",
  "uniprot_sptrembl",
  "rnacentral",
  "refseq_mrna",
  "refseq_mrna_predicted",
  "refseq_ncrna",
  "refseq_ncrna_predicted",    
  "rfam"     
)

OUTPUT_RNA_ATTRIBUTES = "/home/diogo/Documents/RAINET_data/BioMart/RNA_ATTRIBUTES_ensembl82.tsv"
OUTPUT_RNA_XREF_ATTRIBUTES = "/home/diogo/Documents/RAINET_data/BioMart/RNA_XREF_ATTRIBUTES_ensembl82.tsv"

if (file.exists(OUTPUT_RNA_ATTRIBUTES)) {
  warning("The following output file already exists, this script will append data to it: ",OUTPUT_RNA_ATTRIBUTES)
}

BATCH_SIZE = 500 #500 is the advised maximum limit of filter items from BioMart

# Establishing server connection
##########################################################################

ensembl = useMart(host= ENSEMBL_URL, biomart= ENSEMBL_BIOMART, dataset= BIOMART_DATASET)
paste("Accessing host:",ensembl@host)
#listMarts(host="www.ensembl.org")
#listDatasets(ensembl)
#attributes = listAttributes(ensembl)

# Getting all transcript IDs
##########################################################################

allRNAs = getBM(attributes= BIOMART_CENTRAL_ATTRIBUTE, mart = ensembl)
#allRNAs = getBM(attributes= BIOMART_CENTRAL_ATTRIBUTE, filters = "chromosome_name", values = "13", mart = ensembl) #testing set with 3241 transcripts
nRNAs = nrow(allRNAs)
paste("Total number of RNAs:",nRNAs)

# Retrieve features in transcript batches, write to output
##########################################################################
#Note: transcript features have 1-to-1 correspondence, e.g. one transcript only has one associated gene or genomic coordinate.
#Wanted output example (.tsv):
#P62258	1433E_HUMAN	14-3-3 protein epsilon (14-3-3E)	YWHAE		Homo sapiens (Human)	255		PF00244;	SM00101;
#Q04917	1433F_HUMAN	14-3-3 protein eta (Protein AS1)	YWHAH	YWHA1	Homo sapiens (Human)	246		PF00244;	SM00101;

paste("Retrieving transcript models..")

#Write header
header = capture.output(cat(RNA_ATTRIBUTES,sep ="\t") )
write.table(header,file = OUTPUT_RNA_ATTRIBUTES,quote = FALSE,sep="\t", append = TRUE,  row.names = FALSE, col.names = FALSE)

#Approach: for batches of transcripts, query all wanted attributes of that transcript.
j = 0
rnaAttStore = data.frame()
for(i in seq(from = BATCH_SIZE, to = nRNAs+BATCH_SIZE, by = BATCH_SIZE))
{
  rnaBatch <- allRNAs$ensembl_transcript_id[(j+1):min(i,nRNAs)] #parenthesis on (j+1) are very important for R!
  rnaAtt = getBM(attributes=RNA_ATTRIBUTES, filters = BIOMART_CENTRAL_ATTRIBUTE, values = rnaBatch,mart = ensembl)
  rnaAttStore = rbind(rnaAttStore,rnaAtt)
  print(paste("Processed..",round((min(i,nRNAs) *100 / nRNAs),digits = 2 ),"%"))
  j = i
}

write.table(rnaAttStore,file = OUTPUT_RNA_ATTRIBUTES, quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE, append = TRUE)
paste(OUTPUT_RNA_ATTRIBUTES," file written.")

# Retrieve cross references in transcript batches, write to output
##########################################################################
#Note: two constraints in this section: 1) can only retrieve maximum 3 types of external ID on each query; 2) correspondence between ensembl_transcript_id and other xrefs is 1-to-many
#Wanted output example (.tsv):
#P31946	GI	21328448
#P31946	UniRef100	UniRef100_P31946

paste("Retrieving cross references..")

#Approach: for batches of transcripts, query one xref attribute at with each query.
j = 0
rnaAttStore = data.frame()
for(i in seq(from = BATCH_SIZE, to = nRNAs+BATCH_SIZE, by = BATCH_SIZE))
{ ## for each batch of transcripts
  rnaBatch <- allRNAs$ensembl_transcript_id[(j+1):min(i,nRNAs)] 
  
  ## make one query for each xref attribute
  for(k in 1:length(RNA_CROSS_REFERENCE_ATTRIBUTES)) 
  { #for each attribute
    
    #query several transcripts for a single attribute (plus ID)
    rnaAtt = getBM(attributes=c(BIOMART_CENTRAL_ATTRIBUTE,RNA_CROSS_REFERENCE_ATTRIBUTES[k]), filters = BIOMART_CENTRAL_ATTRIBUTE, values = rnaBatch, mart = ensembl)

    #filter out NA and empty lines
    filt = subset(x = rnaAtt, subset = !is.na(rnaAtt[,2]) & rnaAtt[,2] != "") 
    
    #write into file only the items with a value after filtering
    if (nrow(filt)>0){ 
      #adding column with the xref DB name
      filt$db = RNA_CROSS_REFERENCE_ATTRIBUTES[k]
      #swap column
      filt = data.frame(filt[,1],filt[,3],filt[,2])
      rnaAttStore = rbind(rnaAttStore,filt)

    } 
  }

  print(paste("Processed..",round((min(i,nRNAs) *100 / nRNAs),digits = 2 ),"%"))
  j = i  
}

write.table(rnaAttStore, file = OUTPUT_RNA_XREF_ATTRIBUTES, quote = FALSE,sep="\t", row.names = FALSE, col.names = FALSE, append = TRUE)
paste(OUTPUT_RNA_XREF_ATTRIBUTES," file written.")

print ("FINISHED!")

