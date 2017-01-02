
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
metric = "Syn" # change here the wanted metric
dataset1 = datasetContinuous[datasetContinuous$Metric == metric]
table1 = all_vs_all_tests(dataset1, "Value", "Group", verbose = 1)

### Get list of enriched lncRNAs in clusters c1,c2 and c3
clEnrich = datasetDiscrete[datasetDiscrete$Group == "4-Enriched"]
clEnrich[clEnrich$Value == "c1" | clEnrich$Value == "c2" | clEnrich$Value == "c3"]$Gene
#clEnrich[clEnrich$Value == "c1" | clEnrich$Value == "c2" | clEnrich$Value == "c3" | clEnrich$Value == "c4" | clEnrich$Value == "c6"]$Gene



# ##### Network analysis #####
# 
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/Cutoff50/merge/enrichment_network_analysis/topPartners3/lncRNA_group_analysis.tsv"
# 
# dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
# 
# ## Box plot with all numeric metrics at once
# plt1 <- ggplot(dataset, aes(x = Metric, y = Value, fill = Group) )  +
#   geom_boxplot(outlier.shape = NA, position = "dodge") +
#   theme_minimal()
# plt1
# 
# ## All vs All statistics
# metric = "LCneighbours" # change here the wanted metric
# dataset1 = dataset[dataset$Metric == metric]
# table = all_vs_all_tests(dataset1, "Value", "Group", verbose = 1)

