
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

#### Report on data presence ####

inputFile = paste( output_folder, report_rna_expression_data_presence, sep = "/")

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/real/Report/rna_expression_data_presence.tsv"

inputData <- read.table(inputFile, header = TRUE, sep = "\t")

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.7)),
  colhead = list(fg_params=list(cex = 0.8)),
  rowhead = list(fg_params=list(cex = 0.8)))

myt <- gridExtra::tableGrob( inputData, theme = mytheme)

grid.draw(myt)


#### Report on data values ####

#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/rna_expression.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/real/Report/rna_expression.tsv"

inputFile = paste( output_folder, report_rna_expression, sep = "/")

print (paste("INPUT FILE:",inputFile) )

dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

dataLimit = 0.8

print (paste("DISPLAYING DATA WITH XLIM AT PERCENTILE:",dataLimit) )

plt1 <- ggplot(dataset, aes(x = type, y = meanExpression, fill = type ) )  +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() + 
  scale_y_continuous(limits = quantile(dataset$meanExpression, c(0.0, dataLimit))) +
  theme_minimal()
plt1

plt2 <- ggplot(dataset, aes(x = transcriptBiotype, y = meanExpression, fill = type ) )  +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() + 
  scale_y_continuous(limits = quantile(dataset$meanExpression, c(0.0, dataLimit))) +
  theme_minimal()
plt2

grid.arrange( plt1, plt2)


