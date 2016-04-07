
library(ggplot2)
library(gplots)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)
library(data.table)

# NORAD
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/manualLists/NORAD_paper_max_vs_mean.tsv"

dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

datasetMaxOrder <- dataset[order (dataset$catrapid_score, decreasing = TRUE),]
datasetMeanOrder <- dataset[order( dataset$catrapid_score_mean, decreasing = TRUE),] 

fileName = strsplit(inputFile,split="/")[[1]]
fileName = tail(fileName, n=1)

## Scatterplot with the raw points

# correlation, and its significance
correlation = cor(dataset$catrapid_score, dataset$catrapid_score_mean, method = "spearman")
correlationSign = as.numeric(cor.test(dataset$catrapid_score, dataset$catrapid_score_mean, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")

plt2 = ggplot(dataset, aes( x = catrapid_score, y = catrapid_score_mean)) +
  geom_point() +
  geom_smooth(method=lm) + 
  annotate("text", x = max(dataset$catrapid_score), y = max(dataset$catrapid_score_mean), label = correlationText, hjust=1 ) +
  labs(x="max score",y="mean score") +
  ggtitle( fileName)
plt2

#pumilio Q8TB72

