
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)

# function for number of observations 
give.n <- function(x){   return(c(y = -10, label = length(x))) }

# # broad categories files
# inputFile1 = "/home/diogo/testing/RBPDomainScore/lncrna/output/annotated_interactions.tsv"
# inputFile2 = "/home/diogo/testing/RBPDomainScore/mrna/output/annotated_interactions.tsv"

# broad categories files # with count
inputFile1 = "/home/diogo/testing/RBPDomainScore/lncrna/cutoff15/annotated_interactions.tsv"
inputFile2 = "/home/diogo/testing/RBPDomainScore/mrna/cutoff50/annotated_interactions.tsv"

dataset1 <- fread(inputFile1, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2 <- fread(inputFile2, stringsAsFactors = FALSE, header = TRUE, sep="\t")

#adding type before merging datasets
dataset1$type = "lncRNA"
dataset2$type = "mRNA"

head(dataset1)
head(dataset2)

dataset1count = count(dataset1, 'annotation')
dataset2count = count(dataset2, 'annotation')

table = merge( dataset1count, dataset2count, by="annotation", all = TRUE, suffixes = c(".lncRNA", ".mRNA"))
grid.table(table)

#cat(dataset1$uniprotac[dataset1$annotation == "Unclassified"], sep = "\n")
#cat(dataset2$uniprotac[dataset2$annotation == "Unclassified"], sep = "\n")

####################
## Broad categories, Plot distributions
####################

# merge/append the two datasets
mergedDataset = rbind( dataset1, dataset2)

plt1 <- ggplot(dataset1 )  +
  geom_density(data = dataset1, aes(x = count, colour = annotation) ) +
  ggtitle(dataset1$type) +
  theme_minimal()

plt2 <- ggplot(dataset2 )  +
  geom_density(data = dataset2, aes(x = count, colour = annotation) ) +
  ggtitle(dataset2$type) +
  theme_minimal()

grid.arrange(plt1,plt2)

# plt3 <- ggplot(data = mergedDataset, aes(x = annotation, y = count, fill = type))  +
#   geom_boxplot(outlier.shape = NA, position = "dodge" ) +
#   stat_summary(fun.data = give.n, geom = "text", size = 4) +
#   coord_flip() + 
#   theme_minimal()
# plt3

plt3 <- ggplot(data = dataset1, aes(x = annotation, y = count))  +
  geom_boxplot(outlier.shape = NA, position = "dodge", fill = "red" ) +
  stat_summary(fun.data = give.n, geom = "text", size = 4) +
  ggtitle(dataset1$type) +
  coord_flip() + 
  theme_minimal()
plt3

plt4 <- ggplot(data = dataset2, aes(x = annotation, y = count))  +
  geom_boxplot(outlier.shape = NA, position = "dodge", fill = "blue" ) +
  stat_summary(fun.data = give.n, geom = "text", size = 4) +
  ggtitle(dataset2$type) +
  coord_flip() + 
  theme_minimal()
plt4

grid.arrange(plt3,plt4)

####################
## Statistical tests
####################

#test lncRNA classical vs non-classical
lncRNAClassical = dataset1$count[ dataset1$annotation == "Classical"]
lncRNANonClassical = dataset1$count[ dataset1$annotation == "Non-classical"]
lncRNANonRBP = dataset1$count[ dataset1$annotation == "Non-RBP"]
summary( lncRNAClassical)
summary( lncRNANonClassical)
summary( lncRNANonRBP)

ks.test(lncRNAClassical, lncRNANonClassical, alternative = c("two.sided"))
ks.test(lncRNAClassical, lncRNANonRBP, alternative = c("two.sided"))


#test lncRNA classical vs non-classical
mRNAClassical = dataset2$count[ dataset2$annotation == "Classical"]
mRNANonClassical = dataset2$count[ dataset2$annotation == "Non-classical"]
mRNANonRBP = dataset2$count[ dataset2$annotation == "Non-RBP"]

summary( mRNAClassical)
summary( mRNANonClassical)
summary( mRNANonRBP)

ks.test(mRNAClassical, mRNANonClassical, alternative = c("two.sided"))
ks.test(mRNAClassical, mRNANonRBP, alternative = c("two.sided"))

