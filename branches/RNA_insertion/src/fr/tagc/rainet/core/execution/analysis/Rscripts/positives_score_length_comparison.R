
#26-Sep-2016 Script to test idea that interaction signal is diluted when transcript is large. 
#For this reason we pick positive interactions from NPInter/Starbase and compare catRAPID scores with transcript length.

library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

# positives
interactionFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/NPInterPredictionValidation/new_dataset/modified_positive_interactions.txt"
interactionFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/NPInterPredictionValidation/snoRNA_only/modified_methodology_used.tsv"
interactionFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/mRNAs/modified_positive_scores.tsv"
#interactionFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/mRNAs/modified_negative_scores.tsv"
#interactionFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/mRNAs/extraStringent/modified_positive_scores.tsv"

dataset <- fread(interactionFile, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")

# annotation / RNA information
rnaInfoFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/ScoreComparison/RNA_table_RAINET.csv"
rnaInfo <- fread(rnaInfoFile, stringsAsFactors = FALSE, header = TRUE, sep=",")

rnaInfo$transcriptLength = as.numeric( rnaInfo$transcriptLength)

# merge with RNA info
mergedDatasetInfo = merge( dataset, rnaInfo, by="ensembl_id", all.x = TRUE)

#head(mergedDatasetInfo)

# density plot for scores
plt01 <- ggplot(data = mergedDatasetInfo, aes(x = score) )  +
  geom_density( size = 1 ) +
  theme_minimal()
plt01
# density plot for length
plt02 <- ggplot(data = mergedDatasetInfo, aes(x = transcriptLength) )  +
  geom_density( size = 1 ) +
  theme_minimal()
plt02

# to exclude transcripts annotated with having > 1200 nt, which may happen as transcript lengths comes from biomart and not fasta file of catrapid input
mergedDatasetInfo = mergedDatasetInfo[mergedDatasetInfo$transcriptLength <= 1200]


## measure correlation between transcript length and score
correlation = cor(mergedDatasetInfo$transcriptLength, mergedDatasetInfo$score, method = "spearman", use="complete")
correlationSign = as.numeric(cor.test(mergedDatasetInfo$transcriptLength, mergedDatasetInfo$score, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")

plt1 <- ggplot(data = mergedDatasetInfo, aes(x = transcriptLength, y = score) )  +
  geom_point( shape = 1, alpha=1/4 ) +
  geom_smooth( ) +
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
  theme_minimal()
plt1

