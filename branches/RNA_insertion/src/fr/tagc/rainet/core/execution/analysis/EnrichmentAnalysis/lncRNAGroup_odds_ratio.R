
library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/outFile.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/experimental_interactions/outFile.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/rbp_enrichments/outFile.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/rbp_enrichments/outFile_complex_dataset.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/outFile_simple.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/outFile_complex_dataset.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/structure_comparison/outFile.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/structure_comparison/outFile_gencodebasic_background.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# get all different categories  
categories = unique(dataset$ExternalList)

# create a table for each category
for (i in categories){
  grid.newpage()
  grid.table( dataset[dataset$ExternalList == i])
}

# Whole summary
grid.newpage()
grid.table( dataset)
