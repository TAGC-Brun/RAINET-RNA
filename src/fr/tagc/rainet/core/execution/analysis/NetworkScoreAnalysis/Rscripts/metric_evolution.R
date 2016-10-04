
#03-Oct-2016 Script associated to NetworkScoreAnalysis.py, to plot measures of protein closeness with different top partners parameter. 

library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

# positives
resultFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/networkAnalysis/NetworkScoreAnalysis/100tx_produce_plots_extended/combined_results.tsv"


dataset <- fread(resultFile, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")

meltedDataset = melt(dataset, id="TopPartners")

LCNdataset = meltedDataset[ meltedDataset$variable == "LCN_real" |  meltedDataset$variable == "LCN_random"]
SPdataset = meltedDataset[ meltedDataset$variable == "SP_real" |  meltedDataset$variable == "SP_random"]

# LNC differential
dataset$LCN_real_minus_random = dataset$LCN_real - dataset$LCN_random
# SP differential
dataset$SP_real_minus_random = dataset$SP_real - dataset$SP_random


## LCN metric #xlim
plt0 <- ggplot(data = LCNdataset, aes(x = TopPartners, y = value, color = variable) )  +
  geom_point( ) +
  geom_smooth() +
  ylab( "LC neighbours score") +
  xlim( c(2,10)) +
  ylim( c(0,15)) +
  ggtitle( "LC neighbours. Top partners 2 to 10") +
  theme_minimal()
plt0

## LCN metric
plt1.0 <- ggplot(data = LCNdataset, aes(x = TopPartners, y = value, color = variable) )  +
  geom_point( ) +
  geom_smooth() +
  ylab( "LC neighbours score") +
  ggtitle( "LC neighbours. Top partners 2 to 10") +
  theme_minimal()
plt1.0

## LCN metric # all points
plt1.1 <- ggplot(data = dataset, aes(x = TopPartners, y = LCN_real_minus_random) )  +
  geom_point( ) +
  ylab( "LC neighbours score") +
  ggtitle( "Differential of LC neighbours. Top partners 2 to 100") +
  geom_smooth() +
  theme_minimal()
plt1.1

grid.arrange( plt0, plt1.0, plt1.1)

## SP metric # xlim
plt2 <- ggplot(data = SPdataset, aes(x = TopPartners, y = value, color = variable) )  +
  geom_point( ) +
  ylab( "Mean shortest path") +
  xlim( c(2,10)) +
  ggtitle( "Mean shortest path. Top partners 2 to 10") +
  geom_smooth() +
  theme_minimal()
plt2

## SP metric # all points
plt3 <- ggplot(data = SPdataset, aes(x = TopPartners, y = value, color = variable) )  +
  geom_point( ) +
  ylab( "Mean shortest path") +
  ggtitle( "Mean shortest path. Top partners 2 to 100") +
  geom_smooth() +
  theme_minimal()
plt3

## SP metric # all points
plt4 <- ggplot(data = dataset, aes(x = TopPartners, y = SP_real_minus_random) )  +
  geom_point( ) +
  ylab( "Differential shortest path") +
  ggtitle( "Differential of mean shortest path. Top partners 2 to 100") +
  geom_smooth() +
  theme_minimal()
plt4

grid.arrange( plt2, plt3, plt4)



