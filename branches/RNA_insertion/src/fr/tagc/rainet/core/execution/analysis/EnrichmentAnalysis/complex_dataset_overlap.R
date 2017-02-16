
# 15-Feb-2017 Diogo Ribeiro
# Script to plot intra and inter complex-dataset overlap

library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
#source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

#### intra dataset overlap analysis ####

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/complex_dataset_overlap/intra_dataset_results.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/complex_dataset_overlap/interacting_proteins_filter/intra_dataset_results.tsv"

data <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

### data for each dataset
datalist = list()
counter = 0
for (i in unique(data$dataset) ){
  filtData = data[ data$dataset == i]
  
  ## percentage of complexes with at least one high overlap # Bar plot with mean overlap
  highPerc = round( nrow( filtData[filtData$high_annot_overlap > 0]) * 100 / nrow( filtData), digits = 2 )
  ## mean of mean overlap
  meanOver = round(mean( filtData$mean_overlap), digits = 2)

  modI = gsub("\\|", "_vs_", i)
  counter = counter + 1
  datalist[[counter]] <- c(modI, highPerc, meanOver)
  
  plt1 <- ggplot(filtData, aes(x = mean_overlap) )  +
  geom_histogram( binwidth = 0.05) +
  ggtitle(i) +
  theme_minimal()
  print(plt1)
}

results = do.call(rbind, datalist)
colnames(results) <- c("Dataset", "Perc high overlap", "Mean overlap")

# print results as a table
grid.newpage()
grid.table( results)


#### inter dataset overlap analysis ####

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/complex_dataset_overlap/inter_dataset_results.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/complex_dataset_overlap/interacting_proteins_filter/inter_dataset_results.tsv"

data <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

### data for each dataset
datalist = list()
counter = 0
for (i in unique(data$dataset_comparison) ){
  filtData = data[ data$dataset_comparison == i]
  ## percentage of complexes with at least one high overlap # Bar plot with mean overlap
  highPerc = round( nrow( filtData[filtData$high_annot_overlap > 0]) * 100 / nrow( filtData), digits = 2 )
  ## mean of mean overlap
  meanOver = round(mean( filtData$mean_overlap), digits = 2)

  modI = gsub("\\|", "_vs_", i)
  counter = counter + 1
  datalist[[counter]] <- c(modI,highPerc)
  
  #   plt1 <- ggplot( filtData, aes(x = mean_overlap) )  +
  #   geom_histogram( binwidth = 0.05) +
  #   ggtitle(i) +
  #   theme_minimal()
  #   print(plt1)
}

results = do.call(rbind, datalist)
colnames(results) <- c("Dataset comparison", "Perc high overlap", "Mean overlap")

# print results as a table
grid.newpage()
grid.table( results)


