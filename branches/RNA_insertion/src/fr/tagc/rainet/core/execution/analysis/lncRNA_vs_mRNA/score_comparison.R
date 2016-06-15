
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

# catrapid data for each rna
lncRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/output/rnaInteractions.tsv"
mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/mrna/output/rnaInteractions.tsv"

# annotation / RNA information
rnaInfoFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/ScoreComparison/RNA_table_RAINET.csv"
rnaInfo <- fread(rnaInfoFile, stringsAsFactors = FALSE, header = TRUE, sep=",")

dataset1 <- fread(lncRNAData, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")
dataset2 <- fread(mRNAData, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")

# merged lncRNAs and mRNA catrapid data
mergedDataset = rbind( dataset1, dataset2)
stopifnot(nrow(dataset1) + nrow(dataset2) == nrow(mergedDataset))

# merge with RNA info
mergedDatasetInfo = merge( mergedDataset, rnaInfo, by="ensembl_id", all.x = TRUE)

### Plot comparing scores between lncRNA and mRNA
# density plot by category
plt01 <- ggplot(data = mergedDatasetInfo, aes(x = mean_score, colour = type) )  +
  geom_density( size = 1 ) +
#  geom_boxplot( aes(y = mean_score, x = type)) +
  theme_minimal()
# box plot
plt02 <- ggplot(data = mergedDatasetInfo, aes(y = mean_score, x = type) )  +
  geom_boxplot( ) +
  coord_flip() +
  theme_minimal()
grid.arrange(plt01,plt02)
grid.newpage()
table1 = all_vs_all_tests(mergedDatasetInfo, "mean_score", "type", verbose = 1)


## Analyse difference of standard deviation of scores (of transcripts, measure from readCATrapid)

plt5 <- ggplot(data = mergedDatasetInfo, aes(x = std_score, colour = type) )  +
  geom_density( size = 1 ) +
  theme_minimal()
plt6 <- ggplot(data = mergedDatasetInfo, aes(x = type, y = std_score) )  +
  geom_boxplot( ) +
  coord_flip() +
  theme_minimal()
grid.arrange(plt5,plt6)
grid.newpage()
table3 = all_vs_all_tests(mergedDatasetInfo, "std_score", "type", verbose = 1)

#mergedDatasetInfo[ order(mergedDatasetInfo$std_score), ]


### Plot comparing scores between transcript biotypes

# exclude small categories from analysis
MINIMUM_CATEGORY_SIZE = 800
mergedDatasetInfoLess = mergedDatasetInfo
for (i in unique(mergedDatasetInfo$transcriptBiotype) ){
  c = count(mergedDatasetInfoLess$transcriptBiotype == i)$freq[[2]]
  if (c < MINIMUM_CATEGORY_SIZE){    mergedDatasetInfoLess = mergedDatasetInfoLess[mergedDatasetInfoLess$transcriptBiotype != i]  }
}

plt7 <- ggplot(data = mergedDatasetInfoLess, aes(x = transcriptBiotype, y = mean_score) )  +
  geom_boxplot( ) +
  stat_summary(fun.data = give.n, geom = "text", size = 4) +
  coord_flip() +
  theme_minimal()
plt7
grid.newpage()
all_vs_all_tests( mergedDatasetInfoLess, 'mean_score', 'transcriptBiotype', verbose = 1)


# ##### Length control #####
# 
# # ## length data
# # lncRNALengths = "/home/diogo/Documents/RAINET_data/catRAPID/catRAPID_libraries/rna/old/lncrna/seq_lengths/lncrna_lengths.txt"
# # mRNALengths = "/home/diogo/Documents/RAINET_data/catRAPID/catRAPID_libraries/rna/old/cdna/seq_lengths/mrnas_lengths.txt"
# # 
# # dataset1Lengths <- fread(lncRNALengths, stringsAsFactors = FALSE, header = TRUE, sep=" ")
# # dataset2Lengths <- fread(mRNALengths, stringsAsFactors = FALSE, header = TRUE, sep=" ")
# # 
# # # merge lists which contain lengths
# # mergedDatasetLen = rbind(dataset1Lengths, dataset2Lengths)
# # names(mergedDatasetLen) = c("length","ensembl_id")
# # 
# # # add lengths data to dataset
# # mergedDatasetLengths = merge( mergedDataset, mergedDatasetLen, by="ensembl_id", all.x = TRUE)
# 
# mergedDatasetInfo$transcriptLength = as.numeric(mergedDatasetInfo$transcriptLength)
# 
# ## measure correlation between transcript length and score
# correlation = cor(mergedDatasetInfo$transcriptLength, mergedDatasetInfo$mean_score, method = "spearman")
# correlationSign = as.numeric(cor.test(mergedDatasetInfo$transcriptLength, mergedDatasetInfo$mean_score, method = "spearman")$p.value)
# correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")
# 
# plt1 <- ggplot(data = mergedDatasetInfo, aes(x = transcriptLength, y = mean_score) )  +
#   geom_point( shape = 1, alpha=1/4 ) +
#   geom_smooth( ) +
#   annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
#   theme_minimal()
# plt1
# 
# ## filtering so that length distributions match
# 
# # read data
# #rna filter, provide custom lists of transcripts to use
# lncRNAFilt = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/SubsetSimilarLength/list_lncRNAs.txt"
# mRNAFilt = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/SubsetSimilarLength/list_mRNAs.txt"
# 
# dataset1Filt <- fread(lncRNAFilt, stringsAsFactors = FALSE, header = TRUE, sep="\t")
# dataset2Filt <- fread(mRNAFilt, stringsAsFactors = FALSE, header = TRUE, sep="\t")
# 
# # merge lists for filtering
# mergedDatasetFilt = rbind(dataset1Filt, dataset2Filt)
# 
# # Filter previous datasets for the wanted transcripts
# mergedDatasetFiltered = mergedDatasetInfo[ mergedDatasetInfo$ensembl_id %in% mergedDatasetFilt$ensembl_id]
# 
# # validate filter
# nrow(mergedDatasetFiltered) == nrow(mergedDatasetFilt)
# 
# # compare length distribution between lncRNAs and mRNAs
# plt1 <- ggplot(data = mergedDatasetInfo, aes(x = transcriptLength, colour = type) )  +
#   geom_density( size = 2 ) +
#   theme_minimal()
# plt2 <- ggplot(data = mergedDatasetFiltered, aes(x = transcriptLength, colour = type) )  +
#   geom_density( size = 2 ) +
#   theme_minimal()
# grid.arrange(plt1,plt2)
# 
# ## compare scores controlling for length
# plt3 <- ggplot(data = mergedDatasetFiltered, aes(x = mean_score, colour = type) )  +
#   geom_density( size = 1 ) +
#   theme_minimal()
# plt4 <- ggplot(data = mergedDatasetFiltered, aes(x = type, y = mean_score) )  +
#   geom_boxplot( ) +
#   coord_flip() +
#   theme_minimal()
# grid.arrange(plt3,plt4)
# 
# grid.newpage()
# table2 = all_vs_all_tests(mergedDatasetFiltered, "mean_score", "type", verbose = 1)
