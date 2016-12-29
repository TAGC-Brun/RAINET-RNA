
# 29-Dec-2016 Script to plot results of Mukherjee2016 dataset

library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/lncRNA_groups/transcript_to_gene/lncRNA_group_analysis.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

datasetSyn = dataset[dataset$Metric == "Syn"]

datasetSyn


plt1 <- ggplot(dataset, aes(x = Metric, y = Value, fill = Group) )  +
  geom_boxplot() +
  theme_minimal()
plt1

