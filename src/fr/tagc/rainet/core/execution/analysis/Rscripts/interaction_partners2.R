
# Script to plot distribution of interaction partners. From *Interactions.tsv files

library(gplots)
library(RColorBrewer)
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

## Read rnaInteractions.tsv
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/cutoff50/lnc_expr1.58_cutoff50/rnaInteractions.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/lnc_expression_filter_1.58/lnc_expr1.58_cutoff15/rnaInteractions.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")
dataset = dataset[dataset$count != "NA"]

## Read proteinInteractions.tsv
inputFile2 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/cutoff50/lnc_expr1.58_cutoff50/proteinInteractions.tsv"

dataset2 <- fread(inputFile2, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")
dataset2 = dataset2[dataset2$count != "NA"]


##### RNA ######

# for cumsum
sortedDataset = dataset[order(dataset$count)]

# for percentage distribution
mean(dataset$count)

#Plot distribution proteins per RNA

#dataset$countPerc = (dataset$count * 100) / 15974 # nrow(dataset2)

plt1 <- ggplot( dataset, aes(x = count)) + 
  geom_histogram( binwidth = 50) + 
  ylab("# lncRNAs") +
  # scale_y_log10() +
#   ylim(c(0,110)) +
  xlab("# Protein (targets)") +
  ggtitle("Distribution of interactions per lncRNA") +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=4, label = paste("     Median:",median( dataset$count)) ) +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=6, label = paste("Mean:",round( mean( dataset$count), 2 ))) +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=8, label = paste("       Std:",round( sd( dataset$count), 2) )) +
  theme_minimal()
plt1

## Plot cumulative sum

plt2 <- ggplot( sortedDataset, aes(x = 1:nrow(dataset), y = cumsum(count))) + 
  geom_line( ) + 
  ylab("Cumulative sum of target count") +
  xlab("lncRNAs sorted by target count") +
  theme_minimal()
plt2

grid.arrange(plt1, plt2)

##### Proteins ######

# for cumsum
sortedDataset2 = dataset2[order(dataset2$count)]

## Plot distribution RNAs per proteins

mean(dataset2$count)

dataset2$countPerc = (dataset2$count) * 100 / 22960 #/ nrow(dataset)

# plt3 <- ggplot( dataset2, aes(x = countPerc)) + 
#   geom_histogram( binwidth = 0.5) + 
#   ylab("Protein frequency") +
#   xlab("% RNA targets") +
#   theme_minimal()
# plt3

plt3 <- ggplot( dataset2, aes(x = count)) + 
  geom_histogram( binwidth = 25) + 
  ylab("# Proteins") +
  xlab("# lncRNA (targets)") +
  ggtitle("Distribution of interactions per protein") +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=4, label = paste("     Median:",median( dataset2$count)) ) +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=6, label = paste("Mean:",round( mean( dataset2$count), 2 ))) +
  annotate("text", x = Inf, y = Inf,hjust=2,vjust=8, label = paste("       Std:",round( sd( dataset2$count), 2) )) +
  theme_minimal()
plt3


## Plot cumulative sum

plt4 <- ggplot( sortedDataset2, aes(x = 1:nrow(dataset2), y = cumsum(count))) + 
  geom_line( ) + 
  ylab("Cumulative sum of target count") +
  xlab("Proteins sorted by target count") +
  theme_minimal()
plt4

grid.arrange(plt3, plt4)


nrow(dataset)
nrow(dataset2)
