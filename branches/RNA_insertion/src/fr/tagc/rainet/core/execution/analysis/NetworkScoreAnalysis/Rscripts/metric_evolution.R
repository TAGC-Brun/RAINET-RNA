
#03-Oct-2016 Script associated to NetworkScoreAnalysis.py, to plot measures of protein closeness with different top partners parameter. 

library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

# positives
resultFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/networkAnalysis/NetworkScoreAnalysis/100tx_produce_plots/example_result.txt"
resultFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/networkAnalysis/NetworkScoreAnalysis/100tx_produce_plots/combined_results.tsv"

dataset <- fread(resultFile, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")

meltedDataset = melt(dataset, id="TopPartners")

LCNdataset = meltedDataset[ meltedDataset$variable == "LCN_real" |  meltedDataset$variable == "LCN_random"]
SPdataset = meltedDataset[ meltedDataset$variable == "SP_real" |  meltedDataset$variable == "SP_random"]

## LCN metric
plt1 <- ggplot(data = LCNdataset, aes(x = TopPartners, y = value, color = variable) )  +
  geom_point( ) +
  geom_smooth() +
  theme_minimal()
plt1

## SP metric
plt2 <- ggplot(data = SPdataset, aes(x = TopPartners, y = value, color = variable) )  +
  geom_point( ) +
  geom_smooth() +
  theme_minimal()
plt2

