
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)

inputFile = "/home/diogo/testing/RBPDomainScore/output/annotated_interactions.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

head(dataset)

counts = count(dataset, 'annotation')
counts

plt1 <- ggplot(dataset, aes(x = annotation, y = score ) )  +
  geom_boxplot() +
  coord_flip() + 
  theme_minimal()
plt1

