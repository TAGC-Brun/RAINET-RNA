
# 20-Jan-2017 Diogo Ribeiro
# Script to plot proportion of proteins with a certain annotation, between groups of proteins based on lncRNA/mRNA binding ratio.

library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
#source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

# Arguments
inputFileOdds = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/mrna_vs_lncrna/t_test/protein_target_ratio_ttest_RBP_only.out"
#inputFileOdds ="/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/mrna_vs_lncrna/RBP_only/protein_target_ratio_rbp_only_cutoff50_cutoff50.tsv"
#inputFileInfo = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_vs_mRNA/mrna_vs_lncrna/t_test/target_TF/annotated_interactions.tsv"
inputFileInfo = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/backup/lncrna/broad_TFs/annotated_interactions.tsv"

topNumber = 50
oddsType = "t_test_statistic" #"fisher_odds" #

# Read input files
datasetOdds <- fread(inputFileOdds, stringsAsFactors = FALSE, header = TRUE, sep="\t")
datasetInfo <- fread(inputFileInfo, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Process input files
datasetSorted = datasetOdds[ order( datasetOdds[[oddsType]])]

# get the 3 datasets, top, bottom and the ones that not significant (no preference for lncRNA or mRNA)
topDataset = head( datasetSorted, n = topNumber)
bottomDataset = tail( datasetSorted, n = topNumber)
middleDataset = datasetSorted[ datasetSorted$significant == 0]
wholeDataset = datasetSorted

# combine all 3 datasets
topDataset$dataset = "top mRNA>lncRNA"
bottomDataset$dataset = "top lncRNA>mRNA"
middleDataset$dataset = "lncRNA = mRNA"
wholeDataset$dataset = "all proteins"
mergeDataset = rbind( topDataset, bottomDataset, middleDataset, wholeDataset)

# add annotation information on dataset
mergeDatasetInfo <- merge( mergeDataset, datasetInfo, by.x=c("proteinID"), by.y=c("uniprotac"))
# mergeDatasetInfo contain duplicates, since datasets are not mutually exclusive

## Bar plot with proportion per categories
plt1 <- ggplot(mergeDatasetInfo, aes(x = dataset, fill = annotation) )  +
  geom_bar( position = "fill") +
  ylab("Proportion of RPBs") + 
  theme_minimal()
plt1


### Fisher exact tests 
category = "Classical" # change here the wanted metric
group = "top lncRNA>mRNA"

categoryDF = mergeDatasetInfo[mergeDatasetInfo$annotation == category]
categoryAll = categoryDF[categoryDF$dataset == "all proteins"]
groupOverlap = categoryDF[categoryDF$dataset == group]
groupNonOverlap = mergeDatasetInfo[mergeDatasetInfo$dataset == group]

q1 = nrow(groupOverlap)
q2 = nrow(groupNonOverlap) - q1
q3 = nrow(categoryDF) - q1
q4 = nrow(mergeDatasetInfo) - q1 - q2 - q3

fisher.test( matrix(c(q1,q2,q3,q4), nrow = 2))


