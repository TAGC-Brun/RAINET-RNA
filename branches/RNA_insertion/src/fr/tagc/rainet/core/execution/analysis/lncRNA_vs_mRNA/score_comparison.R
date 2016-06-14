
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")


# catrapid data for each rna
lncRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/catrapid_parsing/lncrna/output/rnaInteractions.tsv"
mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/catrapid_parsing/mrna/output/rnaInteractions.tsv"

# # catrapid data for each protein
# lncRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/catrapid_parsing/lncrna/output/proteinInteractions.tsv"
# mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/catrapid_parsing/mrna/output/proteinInteractions.tsv"

dataset1 <- fread(lncRNAData, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")
dataset2 <- fread(mRNAData, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")

#adding type to distinguish datasets
dataset1$type = "lncRNA"
dataset2$type = "mRNA"

mergedDataset = rbind( dataset1, dataset2)

# density plot by category
plt01 <- ggplot(data = mergedDataset, aes(x = mean_score, colour = type) )  +
  geom_density( size = 1 ) +
#  geom_boxplot( aes(y = mean_score, x = type)) +
  theme_minimal()
plt01

# box plot
plt02 <- ggplot(data = mergedDataset, aes(y = mean_score, x = type) )  +
  geom_boxplot( ) +
  theme_minimal()
plt02

grid.arrange(plt01,plt02)


##### Length control #####

## length data
lncRNALengths = "/home/diogo/Documents/RAINET_data/catRAPID/catRAPID_libraries/rna/old/lncrna/seq_lengths/lncrna_lengths.txt"
mRNALengths = "/home/diogo/Documents/RAINET_data/catRAPID/catRAPID_libraries/rna/old/cdna/seq_lengths/mrnas_lengths.txt"

dataset1Lengths <- fread(lncRNALengths, stringsAsFactors = FALSE, header = TRUE, sep=" ")
dataset2Lengths <- fread(mRNALengths, stringsAsFactors = FALSE, header = TRUE, sep=" ")

# merge lists which contain lengths
mergedDatasetLen = rbind(dataset1Lengths, dataset2Lengths)
names(mergedDatasetLen) = c("length","ensembl_id")

# add lengths data to dataset
mergedDatasetLengths = merge( mergedDataset, mergedDatasetLen, by="ensembl_id", all = TRUE)


## measure correlation between transcript length and score
correlation = cor(mergedDatasetLengths$length, mergedDatasetLengths$mean_score, method = "spearman")
correlationSign = as.numeric(cor.test(mergedDatasetLengths$length, mergedDatasetLengths$mean_score, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")

plt1 <- ggplot(data = mergedDatasetLengths, aes(x = length, y = mean_score) )  +
  geom_point( shape = 1, alpha=1/4 ) +
  geom_smooth( ) +
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
  theme_minimal()
plt1

grid.newpage()
table1 = all_vs_all_tests(mergedDatasetLengths, "mean_score", "type", verbose = 1)


## filtering so that length distributions match

# read data
#rna filter, provide custom lists of transcripts to use
lncRNAFilt = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/SubsetSimilarLength/list_lncRNAs.txt"
mRNAFilt = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/SubsetSimilarLength/list_mRNAs.txt"

dataset1Filt <- fread(lncRNAFilt, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2Filt <- fread(mRNAFilt, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# merge lists for filtering
mergedDatasetFilt = rbind(dataset1Filt, dataset2Filt)

# Filter previous datasets for the wanted transcripts
mergedDatasetFiltered = mergedDatasetLengths[ mergedDatasetLengths$ensembl_id %in% mergedDatasetFilt$ensembl_id]

# validate filter
nrow(mergedDataset) > nrow(mergedDatasetFiltered)
nrow(mergedDatasetFiltered) == nrow(mergedDatasetFilt)

# compare length distribution between lncRNAs and mRNAs
plt1 <- ggplot(data = mergedDatasetLengths, aes(x = length, colour = type) )  +
  geom_density( size = 2 ) +
  theme_minimal()
plt1
plt2 <- ggplot(data = mergedDatasetFiltered, aes(x = length, colour = type) )  +
  geom_density( size = 2 ) +
  theme_minimal()
plt2

grid.arrange(plt1,plt2)

## compare scores controlling for length

plt3 <- ggplot(data = mergedDatasetFiltered, aes(x = mean_score, colour = type) )  +
  geom_density( size = 1 ) +
  theme_minimal()
plt3

plt4 <- ggplot(data = mergedDatasetFiltered, aes(x = type, y = mean_score) )  +
  geom_boxplot( ) +
  theme_minimal()
plt4

grid.arrange(plt3,plt4)

grid.newpage()
table2 = all_vs_all_tests(mergedDatasetFiltered, "mean_score", "type", verbose = 1)
# grid.newpage()
# table2 = all_vs_all_tests(mergedDatasetFiltered, "std_score", "type", verbose = 1)
# grid.arrange(
#   tableGrob(table1),
#   tableGrob(table2),
#   nrow=2)


### ANALYSE STD
# of transcript,s measure from readCATARAPID

