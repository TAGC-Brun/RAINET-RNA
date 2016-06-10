
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

####################
## Parameters
####################

# low complexity data for each protein
lowComplexityData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/low_complexity/segmasker_run/perc_low_complexity.txt"

# catrapid data for each protein
# lncRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/lncrna/stawiski_TF_cutoff15/annotated_interactions.tsv"
# mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/mrna/stawiski_TF_cutoff50/annotated_interactions.tsv"

# catrapid data for each protein
lncRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/lncrna/stawiski_TF/annotated_interactions.tsv"
mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/mrna/stawiski_TF/annotated_interactions.tsv"

# Choose here which metric to use (name of column to use)
metricToUse1 = "mean_score" # for catrapid file
metricToUse2 = "perc_low" # for low complexity file
annotCol = "annotation" # protein annotation/ category column

# Constants
MINIMUM_CATEGORY_SIZE = 5 # annotations/categories with less numbers than this will not be analysed

####################
## Initialisation
####################

dataset1 <- fread(lncRNAData, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")
dataset2 <- fread(mRNAData, stringsAsFactors = FALSE, header = TRUE, sep="\t", na.strings="NA")
lowcomplexityDataset <- fread(lowComplexityData, stringsAsFactors = FALSE, header = TRUE, sep="\t")

#adding type to distinguish datasets
dataset1$type = "lncRNA"
dataset2$type = "mRNA"

# add low complexity data to score/annotation datasets
dataset1 = merge(x = dataset1, y = lowcomplexityDataset, by = "uniprotac", all.x = TRUE)
dataset2 = merge(x = dataset2, y = lowcomplexityDataset, by = "uniprotac", all.x = TRUE)

# using variable name to get data of interest
dataset1$metricToUse1 = dataset1[[metricToUse1]]
dataset1$metricToUse2 = dataset1[[metricToUse2]]
dataset2$metricToUse1 = dataset2[[metricToUse1]]
dataset2$metricToUse2 = dataset2[[metricToUse2]]

# exclude small categories from analysis
for (i in table$annotation){
  count = (table$freq.lncRNA[table$annotation == i])  
  if (count < MINIMUM_CATEGORY_SIZE){
    dataset1 = dataset1[dataset1$annotation != i]
    dataset2 = dataset2[dataset2$annotation != i]
  }
}

# to confirm if both RNA datasets use same protein dataset
#table1 = merge( count(dataset1, annotCol), count(dataset2, annotCol),all = TRUE, by=annotCol, suffixes = c(".lncRNA", ".mRNA"))

# create table with protein numbers per annotation
table1 = count(dataset1, annotCol)
names(table1) = c(annotCol,"freq")
grid.table(table1)

####################
## Low complexity analysis
####################
# low complexity only 
# Note: only need to use one dataset, proteins are the same between datasets)

# low complexity % histogram
plt1 <- ggplot(data = dataset1, aes_string(x = metricToUse2) )  +
  geom_histogram( size = 2, binwidth = 1 ) +
  ggtitle( "low complexity, whole dataset") +
  theme_minimal()

# density plot by category
plt2 <- ggplot(data = dataset1, aes_string(x = metricToUse2, colour = "annotation") )  +
  geom_density( size = 1 ) +
  ggtitle( "low complexity, categories") +
  theme_minimal()

grid.arrange( plt1, plt2)

# all vs all tests
all_vs_all_tests( dataset1, metricToUse2, annotCol, verbose = 1)



### Proteins with or without low complexity region

# plots for each category/annotation of proteins

filtAnnot = "*"
plot_filt_dataset(dataset1, dataset2, filtAnnot)

filtAnnot = "RBP"
plot_filt_dataset(dataset1, dataset2, filtAnnot)

filtAnnot = "Non-binding"
plot_filt_dataset(dataset1, dataset2, filtAnnot)

filtAnnot = "Stawiski_control"
plot_filt_dataset(dataset1, dataset2, filtAnnot)


### Low terminus analysis

plt1 <- ggplot(data = dataset1, aes(x = metricToUse1, colour = as.factor(low_terminus)) )  +
  geom_density( size = 1 ) +
  theme_minimal()
plt1


# ###  correlation low complexity vs score
# 
# correlation = cor(dataset1$metricToUse1, dataset1$metricToUse2, method = "spearman")
# correlationSign = as.numeric(cor.test(dataset1$metricToUse1, dataset1$metricToUse2, method = "spearman")$p.value)
# correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")
# 
# plt1 <- ggplot(data = dataset1, aes_string(x = metricToUse1, y = metricToUse2, colour = "annotation"))  +
# #  stat_bin2d(bins = 20) +
#   geom_point( shape = 1, alpha=1/4 ) +
#   geom_smooth( ) +
#   ggtitle(dataset1$type) +
#   annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
#   theme_minimal()
# plt1
# 
# # plt1 <- ggplot(data = dataset1, aes_string(x = metricToUse1, y = metricToUse2) )  +
# #   stat_bin2d(bins = 20) +
# #   ggtitle(dataset1$type) +
# #   annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
# #   theme_minimal()
# # plt1
