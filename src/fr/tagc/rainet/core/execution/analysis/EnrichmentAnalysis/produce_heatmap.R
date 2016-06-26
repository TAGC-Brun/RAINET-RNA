
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/parsing/NetworkModuleR1000/enrichment_results_filtered_matrix.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/parsing/KEGGR100/enrichment_results_filtered_matrix.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")

rownames = dataset$RNAs

dataset <- subset( dataset, select = -RNAs )

rownames(dataset) = rownames

mat_data <- data.matrix(dataset)

nrow(mat_data)
ncol(mat_data)

#colMeans(mat_data)
# apply(mat_data, 1, mean)

par(cex.main=.6)# change the font size of title of plot

my_palette <- colorRampPalette(c("#67a9cf", "#f7f7f7"))(n = 299)

heatmap.2(mat_data,
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
#  dendrogram="none",
  main = "Enrichment of lincRNA interactions on network modules",
  margins=c(1,1),
  distfun = dist,
  hclustfun = hclust,
  key =T, # remove color key
  key.title = "",
  key.xlab = "corrected p-value",
  key.ylab = NULL,
  col=my_palette,       # use on color palette defined earlier 
#   colsep=c(1:6),rowsep=(1:62),
#   sepwidth=c(0.05,0.05), sepcolor="white",
#   Colv="NA",    # turn off column clustering
#   Rowv="NA", #turn off row clustering
#          na.color="blue",
#          symbreaks = min(mat_data, na.rm=TRUE),
)


# dist.pear <- function(x) as.dist(1-cor(t(x)))
# hclust.ave <- function(x) hclust(x, method="average")
# heatmap.2(mat_data, trace="none", distfun=dist.pear, hclustfun=hclust.ave)

# 
# heatmap.2(df,
#           #  cellnote = mat_data,  # same data set for cell labels
#           main = "Test", # heat map title
#           #  notecol="black",      # change font color of cell labels to black
#           density.info="none",  # turns off density plot inside color legend
#           trace="none",         # turns off trace lines inside the heat map
#           margins =c(8,13),     # widens margins around plot
#           col=my_palette,       # use on color palette defined earlier 
#           breaks=col_breaks,    # enable color transition at specified limits
#           #  dendrogram="none",
#           dendrogram="row",
#           #  scale="none",
#           Colv="NA",    # turn off column clustering
#           Rowv = rep_tree_d,
#           #  Rowv="NA",
#           # RowSideColors = c(    # grouping row-variables into different
#           #     rep("green", 37),   # categories, Measurement 1-3: green
#           #     rep("black", 54)),
#           RowSideColors = categ,
#           cexRow = 0.6,
#           cexCol = 0.7,
#           offsetCol = 0.0,
#           offsetRow = 0.1,
#           srtCol = 0,
#           adjCol = c(NA,0),
#           #lwid=c(3, 4, 3 ),
# )         
# labRow = NULL,
# labCol = NULL,
# xlab = "Network modules",
# ylab = "lincRNAs",

# source("https://bioconductor.org/biocLite.R")
# biocLite("made4")
# 
# library("made4")
# require(made4)
# heatplot(mat_data)

