
library(ggplot2)
library(gplots)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)
library(data.table)

#### RNA-Protein interaction matrix heatmap ####

#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/interaction_score_matrix.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/real/Report/interaction_score_matrix.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

rownames = dataset$RNAs

dataset <- subset( dataset, select = -RNAs )

rownames(dataset) = rownames

mat_data <- data.matrix(dataset)

nrow(mat_data)
ncol(mat_data)

# colMeans(mat_data)
# apply(mat_data, 1, mean)

heatmap.2(mat_data,
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          )

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


####### Scatterplots of score sum and protein partners and length

#### RNA ####
inputFile = "/home/diogo/testing/matrix_file/tx_score_length.tsv"
datasetRNA <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

correlation = cor(datasetRNA$scoreSum, datasetRNA$txLength, method = "spearman")
correlationSign = as.numeric(cor.test(datasetRNA$scoreSum, datasetRNA$txLength, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")

plt1 = ggplot(datasetRNA, aes( x = txLength, y = scoreSum)) +
  geom_point(shape=1) + 
  geom_smooth(method=lm) + 
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  )

correlation = cor(datasetRNA$scoreSum, datasetRNA$numProtAbove, method = "spearman")
correlationSign = as.numeric(cor.test(datasetRNA$scoreSum, datasetRNA$txLength, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")

plt2 = ggplot(datasetRNA, aes( x = numProtAbove, y = scoreSum)) +
  geom_point(shape=1) + 
  geom_smooth(method=lm) + 
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  )

grid.arrange( plt1, plt2)

#### Protein ####

inputFile = "/home/diogo/testing/matrix_file/prot_score.tsv"

datasetProt <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
head(datasetProt, 10)
tail(datasetProt, 10)

plt3 <- ggplot( datasetProt, aes(x = numRNAAbove)) + 
  geom_histogram(binwidth = 50, fill="black", colour="black") + 
  theme_minimal() +
  xlab( "Number of interacting RNAs" ) +
  ylab( "Protein count")
plt3

plt4 <- ggplot() + 
  geom_density( data = datasetProt, aes(x= scoreSTD), color = "blue") + 
  geom_density( data = datasetRNA, aes(x= scoreSTD), color = "red") + 
  theme_minimal() +
  xlab( "Score Stdev") +
  ylab( "Density")
plt4


