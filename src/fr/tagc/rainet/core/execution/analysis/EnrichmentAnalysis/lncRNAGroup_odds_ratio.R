
library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/outFile.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")


# get all different categories  
categories = unique(dataset$ExternalList)

# create a table for each category
for (i in categories){
  grid.newpage()
  grid.table( dataset[dataset$ExternalList == i])
}


