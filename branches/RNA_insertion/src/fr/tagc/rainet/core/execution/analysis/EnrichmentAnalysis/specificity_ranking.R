
# Script to plot specificity of enrichments

library(gplots)
library(RColorBrewer)
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/non_redundant/HavugimanaR1000Expr1.0/enrichment_specificity_rank.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/Cutoff50/HavugimanaR10000Expr1.0/enrichment_specificity_rank.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")

##################################
# 2D histogram
##################################

plt1 <- ggplot( dataset, aes(x = transcript_enrichments, y = annot_enrichments)) + 
  stat_bin2d( bin = 20) + 
  xlab("# lncRNA enrichments") +
  ylab("# complex enrichments") +
  ggtitle("2D histogram specificity ranking") +
  xlim( c(1,75)) +
  ylim( c(1,150))
plt1


##################################
# histogram of transcript enrichments
##################################

# get one value per transcript
txDataset = unique( data.frame(dataset$transcriptID, dataset$transcript_enrichments))

plt2 <- ggplot( txDataset, aes(x = dataset.transcript_enrichments)) + 
  geom_histogram( binwidth = 1) +
  ggtitle( "Distribution of number of enrichments per lncRNA") +
  xlab( "Number of enrichments") +
  ylab( "Frequency") + 
  scale_x_continuous( breaks = round(seq(0, max(txDataset$dataset.transcript_enrichments), by = 5))) +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=4, label = paste("     Median:",median( txDataset$dataset.transcript_enrichments)) ) +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=6, label = paste("Mean:",round( mean( txDataset$dataset.transcript_enrichments), 2 ))) +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=8, label = paste("           Std:",round( sd( txDataset$dataset.transcript_enrichments), 2) )) 
plt2


##################################
# histogram of complex enrichments
##################################

# get one value per annotation
annotDataset = unique( data.frame(dataset$annotID, dataset$annot_enrichments))

plt3 <- ggplot( annotDataset, aes(x = dataset.annot_enrichments)) + 
  geom_histogram( binwidth = 5) +
  ggtitle( "Distribution of number of enrichments per complex") +
  xlab( "Number of enrichments") +
  ylab( "Frequency") + 
  scale_x_continuous( breaks = round(seq(0, max(annotDataset$dataset.annot_enrichments), by = 25))) +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=4, label = paste("        Median:",median( annotDataset$dataset.annot_enrichments)) ) +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=6, label = paste("Mean:",round( mean( annotDataset$dataset.annot_enrichments), 2 ))) +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=8, label = paste("       Std:",round( sd( annotDataset$dataset.annot_enrichments), 2) ))
plt3

# 