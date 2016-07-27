
library(data.table)
# library(devtools)
# install_github("ggbiplot", "vqv")
require(ggplot2)
require(grid)
require(gridExtra)
library(ggbiplot)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/test/pca/interaction_score_matrix_altered2.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/mondal2015_ji2016/mondal2015_s5_ji2016/interaction_score_matrix.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/mondal2015_ji2016/annotated_interactions.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/cellular/sample_annotated_interactions.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/ribosomal/sample_annotated_interactions.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

####################
## script options
####################

annotCol = "annotation" # name of the column carrying the annotation information
minimum_category_size = 100 # minimum number of items in category for it to be plotted #0 if no filter

####################
## initialisation
####################

# exclude small categories from analysis
datasetLess = dataset
for (i in unique(dataset$annotation) ){
  c = count(dataset$annotation == i)$freq[[2]]
  if (c < minimum_category_size){    datasetLess = dataset1Less[datasetLess$annotation != i]  }
}
dataset = datasetLess

# use first column as dataframe row names
row.names(dataset) = dataset$RNAs
dataset$RNAs = NULL

# backup dataset, remove annotation column
copyDataset = dataset
dataset[[annotCol]] = NULL
nrow(dataset)
ncol(dataset)

# run PCA:
pca_result <- prcomp(dataset, center = TRUE, scale. = TRUE)

# means and std of each PC
#pca_result$center
#pca_result$scale

# make plot with variance weight for axis
std_dev <- pca_result$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
#scree plot
plot(prop_varex, xlab = "Principal Components",
       ylab = "Proportion of Variance Explained",
       type = "b")
#cumulative scree plot
plot(cumsum(prop_varex), xlab = "Principal Component",
       ylab = "Cumulative Proportion of Variance Explained",
       type = "b")

# # print a summary of the PCA results:
# summary(pca_result)

## biplot with names of PCA, can overload plot if too many
# g <- ggbiplot(pca_result, choices = 1:2, obs.scale = 1, var.scale = 1, ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, circle = TRUE, groups = copyDataset[[annotCol]]) +
#   scale_color_discrete(name = '') +
#   theme(legend.direction = 'horizontal', legend.position = 'top')
# g

listPlot = list()
listPlot[[1]] = ggbiplot(pca_result, choices = 1:2, var.axes = 0, var.scale = 1, alpha = 0.5, ellipse = TRUE,  groups = copyDataset[[annotCol]]) +
 scale_color_discrete(name = '') +
 theme(legend.direction = 'horizontal', legend.position = 'top')
listPlot[[2]] = ggbiplot(pca_result, choices = c(1,3), var.axes = 0, var.scale = 1, alpha = 0.5, ellipse = TRUE,  groups = copyDataset[[annotCol]]) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
listPlot[[3]] = ggbiplot(pca_result, choices = 2:3, var.axes = 0, var.scale = 1, alpha = 0.5, ellipse = TRUE,  groups = copyDataset[[annotCol]]) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
listPlot[[4]] = ggbiplot(pca_result, choices = c(3,4), var.axes = 0, var.scale = 1, alpha = 0.5, ellipse = TRUE,  groups = copyDataset[[annotCol]]) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

grid_arrange_shared_legend(listPlot)

# ### heatmap ###
# 
# library(gplots)
# 
# mat_data <- data.matrix(dataset)
# 
# my_palette <- colorRampPalette(c("#67a9cf", "#f7f7f7"))(n = 20)
# 
# heatmap.2(mat_data,
#           density.info="none",  # turns off density plot inside color legend
#           trace="none",         # turns off trace lines inside the heat map
#           dendrogram="none",
#           col=my_palette,       # use on color palette defined earlier 
# )


# # from https://gist.github.com/thigm85/7689508
# theta <- seq(0,2*pi,length.out = 100)
# circle <- data.frame(x = cos(theta), y = sin(theta))
# plt2 <- ggplot(circle,aes(x,y)) + geom_path()
# loadings <- data.frame(pca_result$rotation,.names = row.names(pca_result$rotation))
# plt2 + geom_text(data=loadings, mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) +
#   coord_fixed(ratio=1) + labs(x = "PC1", y = "PC2")

# library( ade4)
# acp = dudi.pca(dataset)
# 
# cat("\nColumn eigenvalues = ", acp$eig)
# pve <- 100 * acp$eig/sum(acp$eig)
# cat("\nColumn eigenvalues ratios = ", pve)
# cat("\n# of independant variables = ", acp$rank)
# cat("\n# factors kept = ", acp$nf)
# cat("\nCoordinates of variables:\n")
# print(acp$c1)
# cat("\nCoordinates of factors (head):\n")
# print(head(acp$l1))
