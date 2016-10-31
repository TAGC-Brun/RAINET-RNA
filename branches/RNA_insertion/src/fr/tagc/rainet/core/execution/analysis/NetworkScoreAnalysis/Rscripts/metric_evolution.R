
#03-Oct-2016 Script associated to NetworkScoreAnalysis.py, to plot measures of protein closeness with different top partners parameter. 

library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

#
#resultFolder = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/networkAnalysis/NetworkScoreAnalysis/100tx_100R_negative"
resultFolder = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/networkAnalysis/NetworkScoreAnalysis/100tx_shuf_negative"

resultFile = paste( resultFolder, "/combined_results.tsv", sep = "")

dataset <- fread(resultFile, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")

meltedDataset = melt(dataset, id="TopPartners")

LCNdataset = meltedDataset[ meltedDataset$variable == "LCN_real" |  meltedDataset$variable == "LCN_random"]
SPdataset = meltedDataset[ meltedDataset$variable == "SP_real" |  meltedDataset$variable == "SP_random"]

# LNC ratio
dataset$LCN_real_minus_random = dataset$LCN_real / dataset$LCN_random
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
  ggtitle( "LC neighbours") +
  theme_minimal()
plt1.0

## LCN metric # all points
plt1.1 <- ggplot(data = dataset, aes(x = TopPartners, y = LCN_real_minus_random) )  +
  geom_point( ) +
  ylab( "LC neighbours score") +
  ggtitle( "Ratio (real/random) of LC neighbours") +
  geom_smooth() +
  theme_minimal()
plt1.1


## Significant based on LCN metric # all points
plt1.2 <- ggplot(data = dataset, aes(x = TopPartners, y = signLCN) )  +
  geom_bar( stat = "identity" ) +
  ylab( "% significant") +
  ylim( c(0, max(dataset$n_transcripts))) +
  ggtitle( "% significant transcripts per TopPartners") +
  theme_minimal()
plt1.2

grid.arrange( plt0, plt1.0, plt1.1)

grid.arrange( plt1.0, plt1.1)#, plt1.2)


#### SP METRIC

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
  ggtitle( "Mean shortest path. Top partners 2 to Max") +
  geom_smooth() +
  theme_minimal()
plt3

## SP metric # all points
plt4 <- ggplot(data = dataset, aes(x = TopPartners, y = SP_real_minus_random) )  +
  geom_point( ) +
  ylab( "Differential shortest path") +
  ggtitle( "Differential (real - random) of mean shortest path. Top partners 2 to Max") +
  geom_smooth() +
  theme_minimal()
plt4

## Significant based on SP metric # all points
plt5 <- ggplot(data = dataset, aes(x = TopPartners, y = signSP) )  +
  geom_bar( stat = "identity" ) +
  ylab( "# significant") +
  ggtitle( "Number of significant transcripts per TopPartners") +
  theme_minimal()
plt5


grid.arrange( plt2, plt3, plt4)

grid.arrange(  plt3, plt4)#, plt5)

# ### test significance between means of real vs random
# 
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/networkAnalysis/NetworkScoreAnalysis/100tx_1000R/topPartners3/metrics_per_rna.tsv"
# 
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/networkAnalysis/NetworkScoreAnalysis/1000tx_produce_plots/topPartners10/metrics_per_rna.tsv"
# 
# dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")
# 
# mean( dataset$LCneighbours)
# mean( dataset$LCneighboursRandom)
# 
# wilcox_test = wilcox.test( dataset$LCneighbours, dataset$LCneighboursRandom, alternative="greater")
# wilcox_test
# 
# wilcox_test = wilcox.test( dataset$ShortestPath, dataset$ShortestPathRandom, alternative="less")
# wilcox_test
# 

### Combine example transcripts ###
# Evolution of metric for a couple of example transcripts instead of mean

inputFile = paste( resultFolder, "/combined_examples.tsv", sep = "")

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")

meltedDataset = melt(dataset, id="topPartners")

## SP metric # all points
plt10 <- ggplot(data = meltedDataset, aes(x = topPartners, y = value, color = variable) )  +
  geom_line( ) +
  theme_minimal()
plt10


