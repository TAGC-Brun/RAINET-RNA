
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)

# function for number of observations 
give.n <- function(x){   return(c(y = -10, label = length(x))) }

# # domain analysis files
# inputFile1 = "/home/diogo/testing/RBPDomainScore/lncrna/domains/annotated_interactions.tsv"
# inputFile2 = "/home/diogo/testing/RBPDomainScore/mrna/domains/annotated_interactions.tsv"

# domain analysis files + cutoff + TF
inputFile1 = "/home/diogo/testing/RBPDomainScore/lncrna/domains_TF/annotated_interactions.tsv"
inputFile2 = "/home/diogo/testing/RBPDomainScore/mrna/domains_TF/annotated_interactions.tsv"


dataset1 <- fread(inputFile1, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2 <- fread(inputFile2, stringsAsFactors = FALSE, header = TRUE, sep="\t")

#adding type before merging datasets
dataset1$type = "lncRNA"
dataset2$type = "mRNA"

head(dataset1)
head(dataset2)

# dataset1Count = count(dataset1, 'annotation')
# dataset2Count = count(dataset2, 'annotation')
# 
# table = merge( dataset1Count, dataset2Count, by="annotation", all = TRUE, suffixes = c(".lncRNA", ".mRNA"))
# table

####################
## Domain analysis / filtering
####################

dataset1count = count(dataset1, 'annotation')
dataset2count = count(dataset2, 'annotation')

## normalising the two categories
#diff = median( dataset2$mean_score) - median( dataset1$mean_score)
#diff = mean( dataset2$mean_score) - mean( dataset1$mean_score)
#dataset1$mean_score = dataset1$mean_score + diff


## Identify top domains manually
dataset1count[order(dataset1count$freq),]
# list of domains we want to plot (top 15 domains)
wantedDomains = c("RRM_1,RRM_5", "Linker_histone", "Pro_isomerase", "RRM_1,RBM1CTR", "zf-CCCH", "APOBEC_N", "KH_1", "MH2,MH1", "Pkinase", "Helicase_C,DEAD", "RnaseA", "LSM", "WD40", "RRM_1", "Non-binding","TF")

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

plt3 <- ggplot(data = mergedDataset, aes(x = annotation, y = mean_score, fill = type))  +
  geom_boxplot(outlier.shape = NA, position = "dodge" ) +
  stat_summary(fun.data = give.n, geom = "text", size = 4) +
  coord_flip() + 
  theme_minimal()
plt3

####################
## Statistical tests
####################

# #test lncRNA domains
# lncRNARRM_1 = dataset1TopDomains$mean_score[ dataset1TopDomains$annotation == "RRM_1"]
# lncRNALSM = dataset1TopDomains$mean_score[ dataset1TopDomains$annotation == "LSM"]
# lncRNANonRBP = dataset1TopDomains$mean_score[ dataset1TopDomains$annotation == "Non-binding"]
# summary( lncRNARRM_1)
# summary( lncRNALSM)
# summary( lncRNANonRBP)
# 
# ks.test(lncRNARRM_1, lncRNANonRBP, alternative = c("two.sided"))
# ks.test(lncRNALSM, lncRNANonRBP, alternative = c("two.sided"))
# 
# #test mRNA domains
# mRNARRM_1 = dataset2TopDomains$mean_score[ dataset2TopDomains$annotation == "RRM_1"]
# mRNALSM = dataset2TopDomains$mean_score[ dataset2TopDomains$annotation == "LSM"]
# mRNANonRBP = dataset2TopDomains$mean_score[ dataset2TopDomains$annotation == "Non-binding"]
# summary( mRNARRM_1)
# summary( mRNALSM)
# summary( mRNANonRBP)
# 
# ks.test(mRNARRM_1, mRNANonRBP, alternative = c("two.sided"))
# ks.test(mRNALSM, mRNANonRBP, alternative = c("two.sided"))
# 
# # this test is for normalised set
# lncRNAMH2MH1 = dataset1TopDomains$mean_score[ dataset1TopDomains$annotation == "MH2,MH1"]
# mRNAMH2MH1 = dataset2TopDomains$mean_score[ dataset2TopDomains$annotation == "MH2,MH1"]
# ks.test(lncRNAMH2MH1, mRNAMH2MH1, alternative = c("two.sided"))

#### Perform all vs all test ####

categories = unique(dataset1$annotation)

for (i in categories){
  for (j in categories)
    if (i != j){
      set1 = dataset1$metricToUse[ dataset1$annotation == i] 
      set2 = dataset1$metricToUse[ dataset1$annotation == j] 
      print(paste(i, " ", j))
      print(paste("pval:", ks.test(set1, set2, alternative = c("two.sided"))$p.value) )
    }
}


####################
## Number of RNAs above (count)
####################

# domain analysis files + cutoff + TF
inputFile1 = "/home/diogo/testing/RBPDomainScore/lncrna/domains_TF_cutoff15/annotated_interactions.tsv"
inputFile2 = "/home/diogo/testing/RBPDomainScore/mrna/domains_TF_cutoff50/annotated_interactions.tsv"

dataset1 <- fread(inputFile1, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2 <- fread(inputFile2, stringsAsFactors = FALSE, header = TRUE, sep="\t")

dataset1TopDomains = dataset1[ dataset1$annotation %in% wantedDomains]
dataset2TopDomains = dataset2[ dataset2$annotation %in% wantedDomains]
dataset1TopDomains$type = "lncRNA"
dataset2TopDomains$type = "mRNA"

plt3 <- ggplot(data = dataset1TopDomains, aes(x = annotation, y = count))  +
  geom_boxplot(outlier.shape = NA, position = "dodge", fill = "red" ) +
  stat_summary(fun.data = give.n, geom = "text", size = 4) +
  ggtitle(dataset1TopDomains$type) +
  coord_flip() + 
  theme_minimal()

plt4 <- ggplot(data = dataset2TopDomains, aes(x = annotation, y = count))  +
  geom_boxplot(outlier.shape = NA, position = "dodge", fill = "blue" ) +
  stat_summary(fun.data = give.n, geom = "text", size = 4) +
  coord_flip() +
  ggtitle(dataset2TopDomains$type) +
  theme_minimal()

grid.arrange(plt3,plt4)

mRNALSM_count = dataset2TopDomains$count[ dataset2TopDomains$annotation == "LSM"]
mRNANonRBP_count = dataset2TopDomains$count[ dataset2TopDomains$annotation == "Non-binding"]
ks.test(mRNALSM_count, mRNANonRBP_count, alternative = c("two.sided"))

