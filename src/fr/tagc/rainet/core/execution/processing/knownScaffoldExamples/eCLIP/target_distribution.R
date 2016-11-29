
library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

inputFile1 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/eCLIP_interactions/eCLIP_original_interactions/proteinInteractions.tsv"
inputFile2 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/eCLIP_interactions/eCLIP_original_interactions/rnaInteractions.tsv"

dataset1 <- fread(inputFile1, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2 <- fread(inputFile2, stringsAsFactors = FALSE, header = TRUE, sep="\t")

dataset1[order(dataset1$count)]

dataset1$type = "protein"
dataset2$type = "rna"

mean(dataset1$count)


dataset1$uniprotac = NULL
dataset2$ensembl_id = NULL

merged_dataset = rbind(dataset1, dataset2)

plt1 <- ggplot(merged_dataset, aes( x = type, y = count) )  +
  geom_boxplot() +
  coord_flip() +
  theme_minimal()
plt1

