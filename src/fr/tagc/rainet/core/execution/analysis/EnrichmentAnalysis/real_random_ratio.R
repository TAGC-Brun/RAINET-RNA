
library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/CorumR1000Expr/enrichment_per_rna.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/CorumR10000/enrichment_per_rna.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

#head( dataset$n_sign_tests_no_warning, 5)
#head( dataset$avg_n_sign_random_no_warning,5)

dataset$ratio = dataset$n_sign_tests_no_warning / dataset$avg_n_sign_random_no_warning

plt1 <- ggplot(dataset )  +
  geom_density(data = dataset, aes(x = ratio) ) +
  theme_minimal()
plt1

plt1 <- ggplot(dataset )  +
  geom_histogram(data = dataset, aes(x = ratio), binwidth = 0.1 ) +
  xlim(c(0,5)) +
  theme_minimal()
plt1

#tail( dataset[order(dataset$ratio),], 100) 

datasetSign = dataset[dataset$significant == "1"]
nrow(datasetSign)

nrow( datasetSign[datasetSign$ratio >= 2])
nrow( datasetSign[datasetSign$ratio < 2])


