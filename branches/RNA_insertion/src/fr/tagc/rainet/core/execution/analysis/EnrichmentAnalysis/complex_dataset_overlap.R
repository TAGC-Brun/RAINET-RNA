
# 15-Feb-2017 Diogo Ribeiro
# Script to plot intra and inter complex-dataset overlap

library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
#source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/execution/analysis/test_output/complex_dataset_overlap/intra_dataset_results.tsv"

data <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

datasetTag = "BioplexCluster" #WanCluster

# Bar plot with mean overlap
plt1 <- ggplot(data[ data$dataset == datasetTag], aes(x = mean_overlap) )  +
  geom_histogram( binwidth = 0.05) +
  ggtitle(datasetTag) +
  theme_minimal()
plt1

# Bar plot with any overlap
plt2 <- ggplot(data[ data$dataset == datasetTag], aes(x = all_annot_overlap) )  +
  geom_histogram( binwidth = 1) +
  ggtitle(datasetTag) +
  theme_minimal()
plt2

# Bar plot with high overlap
plt3 <- ggplot(data[ data$dataset == datasetTag], aes(x = high_annot_overlap) )  +
  geom_histogram( binwidth = 1) +
  ggtitle(datasetTag) +
  theme_minimal()
plt3

