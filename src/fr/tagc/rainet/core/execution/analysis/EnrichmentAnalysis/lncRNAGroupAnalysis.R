
# 29-Dec-2016 Script to plot results of Mukherjee2016 dataset

library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/lncRNA_groups/transcript_to_gene/lncRNA_group_analysis.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

datasetDiscrete = dataset[dataset$Metric == "Cluster"]
datasetContinuous = dataset[dataset$Metric != "Cluster"]

nrow(datasetDiscrete)
nrow(datasetContinuous)

datasetContinuous$Value <- as.numeric(as.character(datasetContinuous$Value))

## Box plot with all numeric metrics at once
plt1 <- ggplot(datasetContinuous, aes(x = Metric, y = Value, fill = Group) )  +
  geom_boxplot(outlier.shape = NA, position = "dodge") +
  ylim(-2.5,2.5) +
  theme_minimal()
plt1

## Bar plot with cluster categories
plt2 <- ggplot(datasetDiscrete, aes(x = Group, fill = Value) )  +
  geom_bar( position = "fill") +
  theme_minimal()
plt2


### All vs All statistics

metric = "Deg" # change here the wanted metric
dataset1 = datasetContinuous[datasetContinuous$Metric == metric]
table1 = all_vs_all_tests(dataset1, "Value", "Group", verbose = 1)


