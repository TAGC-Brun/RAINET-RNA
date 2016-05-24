
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)

inputFile1 = "/home/diogo/testing/RBPDomainScore/output/annotated_interactions.tsv"
inputFile2 = "/home/diogo/testing/RBPDomainScore/cdna/output/annotated_interactions.tsv"

dataset1 <- fread(inputFile1, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2 <- fread(inputFile2, stringsAsFactors = FALSE, header = TRUE, sep="\t")

head(dataset1)
head(dataset2)

dataset1count = count(dataset1, 'annotation')
dataset2count = count(dataset2, 'annotation')

table = merge( dataset1count, dataset2count, by="annotation", all = TRUE, suffixes = c(".lncRNA", ".mRNA"))
table

plt1 <- ggplot(dataset1, aes(x = annotation, y = score) )  +
  geom_boxplot(outlier.shape = NA, fill = "red") +
  coord_flip() + 
  theme_minimal()
plt1

plt2 <- ggplot(dataset2, aes(x = annotation, y = score ) )  +
  geom_boxplot(outlier.shape = NA, fill = "blue") +
  coord_flip() + 
  theme_minimal()
plt2

grid.arrange(plt1,plt2)

plt3 <- ggplot()  +
  geom_boxplot(outlier.shape = NA, data = dataset1, aes(x = annotation, y = score, colour = "lncRNA" ) ) +
  geom_boxplot(outlier.shape = NA, data = dataset2, aes(x = annotation, y = score, colour = "mRNA" ) ) +
  scale_color_manual(name = "Dataset", values = c("red","blue")) +
  coord_flip() + 
  theme_minimal()
plt3


#> 105+40+123+615+147+56
#[1] 1086

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
summary( mRNAClassical)
summary( mRNANonClassical)
summary( mRNANonRBP)

mRNAClassical_NonClassical = ks.test(mRNAClassical, mRNANonClassical, alternative = c("two.sided"))
mRNAClassical_NonClassical
# for mRNA non classical are significantly higher than classical!


