
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

# # catrapid data for each protein
lncRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/lncrna/stawiski_TF_cutoff15/annotated_interactions.tsv"
mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/mrna/stawiski_TF_cutoff50/annotated_interactions.tsv"

# # # catrapid data for each protein
# lncRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/lncrna/moon_cutoff15/annotated_interactions.tsv"
# mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/mrna/moon_cutoff50/annotated_interactions.tsv"

# catrapid data for each protein
# lncRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/lncrna/stawiski_TF/annotated_interactions.tsv"
# mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/mrna/stawiski_TF/annotated_interactions.tsv"

# Choose here which metric to use (name of column to use)
metricToUse1 = "count" # for catrapid file (mean_score, median_score, std_score, count)
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

# # exclude small categories from analysis
# for (i in table$annotation){
#   count = (table$freq.lncRNA[table$annotation == i])  
#   if (count < MINIMUM_CATEGORY_SIZE){
#     dataset1 = dataset1[dataset1$annotation != i]
#     dataset2 = dataset2[dataset2$annotation != i]
#   }
# }
# 
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


####################
## Low complexity vs score analysis
####################

### Compare scores of proteins with or without low complexity region
# plots for each category/annotation of proteins

filtAnnot = "*"
plot_filt_dataset(dataset1, dataset2, filtAnnot, metricToUse1)

filtAnnot = "RBP"
plot_filt_dataset(dataset1, dataset2, filtAnnot, metricToUse1)

filtAnnot = "Non-binding"
plot_filt_dataset(dataset1, dataset2, filtAnnot, metricToUse1)

filtAnnot = "Stawiski_control"
plot_filt_dataset(dataset1, dataset2, filtAnnot, metricToUse1)

### Low terminus analysis

# idea of separating lc into 3 categories: no LC, terminus LC, central LC

# plot with difference LC locations
plt1 <- ggplot(data = dataset1, aes(x = count, colour = as.factor(tag)) )  +
  geom_density( size = 1 ) +
  theme_minimal()
plt1

# all vs all tests
grid.newpage()
all_vs_all_tests( dataset1, "count", "tag", verbose = 1)


## control for protein length influencing results
# plt1 <- ggplot(data = dataset1, aes_string(x = "total_length", y = "perc_low"))  +
#   geom_point( shape = 1, alpha=1/4 ) +
#   geom_smooth( ) +
#   ggtitle(dataset1$type) +
#   theme_minimal()
# plt1
# 

## correlation between score and % of LC

dataset1Filt = dataset1[dataset1$annotation == "RBP" | dataset1$annotation == "Non-binding" | dataset1$annotation == "TF"]
dataset1Filt = dataset1[dataset1$annotation == "EMFs"]

plt1 <- ggplot(data = dataset1Filt, aes_string(x = "perc_low", y = "count", colour = "annotation"))  +
  geom_point( shape = 1, alpha=1/4 ) +
  geom_smooth( ) +
  ggtitle(dataset1$type) +
  theme_minimal()
plt1



# # instead of using counts, use percentage of transcripts above cutoff
# dataset1$count_perc = dataset1$count * 100 / 8052
# dataset2$count_perc = dataset2$count * 100 / 31378
# 
# plt1 <- ggplot()  +
#   geom_point(data = dataset1, aes_string(x = "count_perc", y = "mean_score"), shape = 1, alpha=1/4, colour = "blue" ) +
#   geom_point(data = dataset2, aes_string(x = "count_perc", y = "mean_score"), shape = 1, alpha=1/4, colour = "red" ) +
# #  coord_flip() +
#   theme_minimal()
# plt1
# 




###  correlation low complexity vs score
# 
# dataset1 = dataset1[dataset1$metricToUse1 != "NA"]
# dataset1 = dataset1[dataset1$annotation == "RBP"]
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
# #   coord_flip() +
#   theme_minimal()
# plt1

# plt1 <- ggplot(data = dataset1, aes_string(x = metricToUse1, y = metricToUse2) )  +
#   stat_bin2d(bins = 20) +
#   ggtitle(dataset1$type) +
#   annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
#   theme_minimal()
# plt1


# ## Starbase data
# 
# mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/forAndreas/protein_interactions.tsv"
# mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/forAndreas/min_reads_100/protein_interactions.tsv"
# 
# dataset2 = dataset2[dataset2$perc_low != "NA"]
# dataset2
# 
# metricToUse1 = "mean_reads"
# metricToUse2 = "perc_low"
# 
# dataset2$metricToUse1 = dataset2[[metricToUse1]]
# dataset2$metricToUse2 = dataset2[[metricToUse2]]
# 
# correlation = cor(dataset2$metricToUse1, dataset2$metricToUse2, method = "spearman")
# correlationSign = as.numeric(cor.test(dataset2$metricToUse1, dataset2$metricToUse2, method = "spearman")$p.value)
# correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")
# 
# plt1 <- ggplot(data = dataset2, aes_string(x = metricToUse1, y = metricToUse2))  +
#   geom_point( shape = 1, alpha=1/4 ) +
#   #  xlim(c(0,8000)) +
#   geom_smooth( ) +
#   annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
#   theme_minimal()
# plt1
# 
# # plot with difference LC locations
# plt1 <- ggplot(data = dataset2, aes(x = count_interactions, colour = as.factor(tag)) )  +
#   geom_density( size = 1 ) +
#   theme_minimal()
# plt1
# 
# # all vs all tests
# grid.newpage()
# all_vs_all_tests( dataset2, "count_interactions", "tag", verbose = 1)

