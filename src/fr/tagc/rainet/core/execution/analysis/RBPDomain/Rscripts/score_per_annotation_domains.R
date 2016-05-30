
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)

# function for number of observations 
give.n <- function(x){   return(c(y = -10, label = length(x))) }

# domain analysis files
inputFile1 = "/home/diogo/testing/RBPDomainScore/lncrna/domains/annotated_interactions.tsv"
inputFile2 = "/home/diogo/testing/RBPDomainScore/mrna/domains/annotated_interactions.tsv"

dataset1 <- fread(inputFile1, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2 <- fread(inputFile2, stringsAsFactors = FALSE, header = TRUE, sep="\t")

#adding type before merging datasets
dataset1$type = "lncRNA"
dataset2$type = "mRNA"

####################
## Domain analysis / filtering
####################

dataset1count = count(dataset1, 'annotation')
dataset2count = count(dataset2, 'annotation')

## Identify top domains manually
dataset1count[order(dataset1count$freq),]
# list of domains we want to plot (top 15 domains)
wantedDomains = c("RRM_1,RRM_5", "Linker_histone", "Pro_isomerase", "RRM_1,RBM1CTR", "zf-CCCH", "APOBEC_N", "KH_1", "MH2,MH1", "Pkinase", "Helicase_C,DEAD", "RnaseA", "LSM", "WD40", "RRM_1", "Non-RBP")

### Filter dataset for having only those top domains
dataset1TopDomains = dataset1[ dataset1$annotation %in% wantedDomains]
dataset2TopDomains = dataset2[ dataset2$annotation %in% wantedDomains]

dataset1TopDomainsCount = count(dataset1TopDomains, 'annotation')
dataset2TopDomainsCount = count(dataset2TopDomains, 'annotation')

# merge/append the two datasets
mergedDataset = rbind( dataset1TopDomains, dataset2TopDomains)

table = merge( dataset1TopDomainsCount, dataset2TopDomainsCount, by="annotation", all = TRUE, suffixes = c(".lncRNA", ".mRNA"))
table

####################
## Plot distributions
####################

plt1 <- ggplot(dataset1TopDomains )  +
  geom_density(data = dataset1TopDomains, aes(x = score, colour = annotation) ) +
  ggtitle(dataset1TopDomains$type) +
  theme_minimal()

plt2 <- ggplot(dataset2TopDomains )  +
  geom_density(data = dataset2TopDomains, aes(x = score, colour = annotation) ) +
  ggtitle(dataset2TopDomains$type) +
  theme_minimal()

grid.arrange(plt1,plt2)

plt3 <- ggplot(data = mergedDataset, aes(x = annotation, y = score, fill = type))  +
  geom_boxplot(outlier.shape = NA, position = "dodge" ) +
  stat_summary(fun.data = give.n, geom = "text", size = 4) +
  coord_flip() + 
  theme_minimal()
plt3

####################
## Statistical tests
####################

#test lncRNA domains
lncRNARRM_1 = dataset1TopDomains$score[ dataset1TopDomains$annotation == "RRM_1"]
lncRNALSM = dataset1TopDomains$score[ dataset1TopDomains$annotation == "LSM"]
lncRNANonRBP = dataset1TopDomains$score[ dataset1TopDomains$annotation == "Non-RBP"]
summary( lncRNARRM_1)
summary( lncRNALSM)
summary( lncRNANonRBP)

ks.test(lncRNARRM_1, lncRNANonRBP, alternative = c("two.sided"))
ks.test(lncRNALSM, lncRNANonRBP, alternative = c("two.sided"))

#test mRNA domains
mRNARRM_1 = dataset2TopDomains$score[ dataset2TopDomains$annotation == "RRM_1"]
mRNALSM = dataset2TopDomains$score[ dataset2TopDomains$annotation == "LSM"]
mRNANonRBP = dataset2TopDomains$score[ dataset2TopDomains$annotation == "Non-RBP"]
summary( mRNARRM_1)
summary( mRNALSM)
summary( mRNANonRBP)

ks.test(mRNARRM_1, mRNANonRBP, alternative = c("two.sided"))
ks.test(mRNALSM, mRNANonRBP, alternative = c("two.sided"))

