
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)

inputFile1 = "/home/diogo/testing/RBPDomainScore/lncrna/output/annotated_interactions.tsv"
inputFile2 = "/home/diogo/testing/RBPDomainScore/cdna/output/annotated_interactions.tsv"

inputFile1 = "/home/diogo/testing/RBPDomainScore/lncrna/output_domains/annotated_interactions.tsv"
inputFile2 = "/home/diogo/testing/RBPDomainScore/cdna/output_domains/annotated_interactions.tsv"

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
table

# merge/append the two datasets
mergedDataset = rbind( dataset1, dataset2)

#cat(dataset1$uniprotac[dataset1$annotation == "Unclassified"], sep = "\n")
#cat(dataset2$uniprotac[dataset2$annotation == "Unclassified"], sep = "\n")

####################
## Plot distributions
####################

plt1 <- ggplot(dataset1 )  +
  geom_density(data = dataset1, aes(x = score, colour = annotation) ) +
  theme_minimal()
plt1

plt2 <- ggplot(dataset2 )  +
  geom_density(data = dataset2, aes(x = score, colour = annotation) ) +
  theme_minimal()
plt2

grid.arrange(plt1,plt2)

plt3 <- ggplot(data = mergedDataset, aes(x = annotation, y = score, fill = type))  +
  geom_boxplot(outlier.shape = NA, position = "dodge" ) +
  coord_flip() + 
  theme_minimal()
plt3

# plt3 <- ggplot()  +
#   geom_boxplot(outlier.shape = NA, data = dataset1, aes(x = annotation, y = score, colour = "lncRNA" ) ) +
#   geom_boxplot(outlier.shape = NA, data = dataset2, aes(x = annotation, y = score, colour = "mRNA" ) ) +
#   scale_color_manual(name = "Dataset", values = c("red","blue")) +
#   coord_flip() + 
#   theme_minimal()
# plt3


####################
## Statistical tests
####################

#test lncRNA classical vs non-classical
lncRNAClassical = dataset1$score[ dataset1$annotation == "Classical"]
lncRNANonClassical = dataset1$score[ dataset1$annotation == "Non-classical"]
lncRNANonRBP = dataset1$score[ dataset1$annotation == "Non-RBP"]
summary( lncRNAClassical)
summary( lncRNANonClassical)
summary( lncRNANonRBP)

lncRNAClassical_NonClassical = ks.test(lncRNAClassical, lncRNANonClassical, alternative = c("two.sided"))
lncRNAClassical_NonClassical

lncRNAClassical_NonRBP = ks.test(lncRNAClassical, lncRNANonRBP, alternative = c("two.sided"))
lncRNAClassical_NonRBP

#for lncRNA no statistical differences

#test lncRNA classical vs non-classical
mRNAClassical = dataset2$score[ dataset2$annotation == "Classical"]
mRNANonClassical = dataset2$score[ dataset2$annotation == "Non-classical"]
mRNANonRBP = dataset2$score[ dataset2$annotation == "Non-RBP"]

summary( mRNAClassical)
summary( mRNANonClassical)
summary( mRNANonRBP)

mRNAClassical_NonClassical = ks.test(mRNAClassical, mRNANonClassical, alternative = c("two.sided"))
mRNAClassical_NonClassical
mRNAClassical_NonRBP = ks.test(mRNAClassical, mRNANonRBP, alternative = c("two.sided"))
mRNAClassical_NonRBP


