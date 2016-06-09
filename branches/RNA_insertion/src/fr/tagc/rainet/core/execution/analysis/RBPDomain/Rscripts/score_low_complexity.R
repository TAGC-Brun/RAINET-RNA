
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

lowComplexityData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/low_complexity/segmasker_run/perc_low_complexity.txt"
lncRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/lncrna/stawiski_cutoff15/annotated_interactions.tsv"
mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/mrna/stawiski_cutoff50/annotated_interactions.tsv"
# lncRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/lncrna/extraMetrics/annotated_interactions.tsv"
# mRNAData = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/RBPDomainScore/mrna/extraMetrics/annotated_interactions.tsv"

dataset1 <- fread(lncRNAData, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2 <- fread(mRNAData, stringsAsFactors = FALSE, header = TRUE, sep="\t")
lowcomplexityDataset <- fread(lowComplexityData, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Choose here which metric to use (name of column to use)
metricToUse1 = "mean_score"
metricToUse2 = "perc_low"
annotCol = "annotation"

#adding type before merging datasets
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

# # check datasets
# head(dataset1)
# head(dataset2)
# table = merge( count(dataset1, 'annotation'), count(dataset2, 'annotation'), by="annotation", all = TRUE, suffixes = c(".lncRNA", ".mRNA"))
# grid.table(table)



####################
## Low complexity analysis
####################

# low complexity only 
#(only need to use one dataset, proteins are the same between datasets)

# low complexity % histogram
plt1 <- ggplot(data = dataset1, aes_string(x = metricToUse2) )  +
  geom_histogram( size = 2, binwidth = 1 ) +
  ggtitle( metricToUse2) +
  theme_minimal()

# density plot by category
plt2 <- ggplot(data = dataset1, aes_string(x = metricToUse2, colour = "annotation") )  +
  geom_density( size = 1 ) +
  ggtitle( metricToUse2) +
  theme_minimal()

grid.arrange( plt1, plt2)


# all vs all tests
annotationCol = annotCol
metricToUse = metricToUse2

source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")
all_vs_all_tests( dataset1, metricToUse, annotationCol)


# low complexity vs score

correlation = cor(dataset1$metricToUse1, dataset1$metricToUse2, method = "spearman")
correlationSign = as.numeric(cor.test(dataset1$metricToUse1, dataset1$metricToUse2, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")

plt1 <- ggplot(data = dataset1, aes_string(x = metricToUse1, y = metricToUse2, colour = "annotation"))  +
#  stat_bin2d(bins = 20) +
  geom_point( shape = 1, alpha=1/4 ) +
  geom_smooth( ) +
  ggtitle(dataset1$type) +
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
  theme_minimal()
plt1


# plt1 <- ggplot(data = dataset1, aes_string(x = metricToUse1, y = metricToUse2) )  +
#   stat_bin2d(bins = 20) +
#   ggtitle(dataset1$type) +
#   annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
#   theme_minimal()
# plt1

