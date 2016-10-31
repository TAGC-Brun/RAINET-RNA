
# Script to plot specificity of enrichments

library(gplots)
library(RColorBrewer)
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/non_redundant/HavugimanaR1000Expr1.0/enrichment_specificity_rank.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")

plt1 <- ggplot( dataset, aes(x = transcript_enrichments, y = annot_enrichments)) + 
  stat_bin2d( bin = 20) + 
  xlab("# lncRNA enrichments") +
  ylab("# complex enrichments") +
  ggtitle("2D histogram specificity ranking") +
  xlim( c(1,75)) +
  ylim( c(1,150))
plt1

# plt1 <- ggplot( dataset, aes(x = transcript_enrichments, y = annot_enrichments)) + 
#   geom_jitter() 
# plt1

# get one value per transcript
txDataset = unique( data.frame(dataset$transcriptID, dataset$transcript_enrichments))
# get one value per annotation
annotDataset = unique( data.frame(dataset$annotID, dataset$annot_enrichments))

plt2 <- ggplot( txDataset, aes(x = dataset.transcript_enrichments)) + 
  geom_histogram( ) +
  theme_minimal()
plt2

plt3 <- ggplot( annotDataset, aes(x = dataset.annot_enrichments)) + 
  geom_histogram( ) +
  theme_minimal()
plt3


