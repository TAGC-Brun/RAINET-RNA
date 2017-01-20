
library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/mrna_vs_lncrna/protein_target_ratio_cutoff50_cutoff50_forR.tsv"

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/mrna_vs_lncrna/t_test/protein_target_ratio_ttest.out"

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/mrna_vs_lncrna/RBP_only/protein_target_ratio_rbp_only_cutoff50_cutoff50.tsv_forR"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

min(dataset$t_test_statistic)
max(dataset$t_test_statistic)
# 

plt1 <- ggplot( dataset, aes(x = dataset$t_test_statistic)) + 
  geom_histogram() + 
  xlab("t test statistic") +
  ylab("# proteins (whole proteome)") +
  theme_minimal()
plt1

