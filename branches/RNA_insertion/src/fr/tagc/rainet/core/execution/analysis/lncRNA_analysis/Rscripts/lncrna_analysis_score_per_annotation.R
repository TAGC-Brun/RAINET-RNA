
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(plyr)

source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

#inputFile1 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/mondal2015_ji2016/nocutoff_s4/annotated_interactions.tsv"
#inputFile1 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/mondal2015_ji2016/cutoff15_s4/annotated_interactions.tsv"

#inputFile1 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/biotypes/nocutoff/annotated_interactions.tsv"
#inputFile1 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/biotypes/cutoff15/annotated_interactions.tsv"

inputFile1 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/carlevaro2016/cellular_cutoff15/annotated_interactions.tsv"
#inputFile1 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/carlevaro2016/ribosomal_cutoff15/annotated_interactions.tsv"

dataset1 <- fread(inputFile1, stringsAsFactors = FALSE, header = TRUE, sep="\t")

####################
## script options
####################

## Choose here which metric to use (name of column to use)
metricToUse = "count" # for catrapid file (mean_score, median_score, std_score, count)
annotCol = "annotation" # name of the column carrying the annotation information
minimum_category_size = 100 # minimum number of items in category for it to be plotted #0 if no filter

####################
## initialisation
####################

# using variable name to get data of interest
dataset1$metricToUse = dataset1[[metricToUse]]

head(dataset1)

dataset1count = count(dataset1, 'annotation')

# exclude small categories from analysis
dataset1Less = dataset1
for (i in unique(dataset1$annotation) ){
  c = count(dataset1$annotation == i)$freq[[2]]
  if (c < minimum_category_size){    dataset1Less = dataset1Less[dataset1Less$annotation != i]  }
}
dataset1 = dataset1Less

####################
## Broad categories, Plot distributions
####################

plt1 <- ggplot(dataset1 )  +
  geom_density(data = dataset1, aes_string(x = metricToUse, colour = "annotation") ) +
  ggtitle(dataset1$type) +
  theme_minimal()
plt1

####################
## Broad categories, Plot distributions as boxplot/violin plot
####################

plt2 <- ggplot(data = dataset1, aes_string(x = "annotation", y = metricToUse))  +
  geom_violin( ) +
  geom_boxplot(outlier.shape = NA, position = "dodge", width = 0.3 ) +
  stat_summary(fun.data = give.n, geom = "text", size = 4) +
  ggtitle(dataset1$type) +
  coord_flip() + 
  theme_minimal()
plt2

# ####################
# ## Statistical tests
# ####################

#### Perform all vs all test ####

grid.newpage()
all_vs_all_tests( dataset1, metricToUse, annotCol)

