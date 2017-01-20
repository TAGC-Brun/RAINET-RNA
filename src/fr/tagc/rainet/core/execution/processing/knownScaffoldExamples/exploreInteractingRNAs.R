
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)
library(plyr)
library(data.table)


inputFile = "/home/diogo/testing/ExploreInteractingRNAs/outFile.txt"

# #lincRNAs only
#inputFile = "/home/diogo/testing/ExploreInteractingRNAs/outFile3.txt"


dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

dataset1 = dataset[dataset$inSet == 0]
dataset2 = dataset[dataset$inSet == 1]

ksTestWhole = ks.test(dataset1$transcriptLength, dataset2$transcriptLength, alternative = c("two.sided"))

plt0 <- ggplot(dataset, aes(x= transcriptLength, colour = as.factor(inSet))) + 
  geom_density() + 
  theme_minimal() +
  xlim( 0,1500) + 
  xlab( "Transcript length") +
  ylab( "Density") +
  annotate("text",  x=Inf, y = Inf, label = paste("# True: ",nrow(dataset2),"\n# False: ",nrow(dataset)-nrow(dataset2)), vjust=1, hjust=1) +
  annotate("text",  x=Inf, y = 0, label = paste("KS test p-val: ", round(ksTestWhole$p.value,2)), vjust=1, hjust=1)
plt0

plt3 <- ggplot( dataset1, aes( x = transcriptBiotype)) + 
  geom_histogram( binwidth = 1) + 
  coord_flip() +
  theme_minimal()

plt4 <- ggplot( dataset2, aes( x = transcriptBiotype)) + 
  geom_histogram( binwidth = 1) + 
  coord_flip() +
  theme_minimal()

grid.arrange(plt3, plt4)


