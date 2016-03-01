
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)


inFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/rna_expression.tsv"

dataset = read.table( inFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

plt1 <- ggplot(dataset, aes(x = type, y = meanExpression, fill = type ) )  +
  geom_boxplot() +
  coord_flip() + 
  theme_minimal()
plt1

plt2 <- ggplot(dataset, aes(x = transcriptBiotype, y = meanExpression, fill = type ) )  +
  geom_boxplot() + #outlier.shape = NA
  coord_flip() + 
  theme_minimal()
plt2

grid.arrange( plt1, plt2)
