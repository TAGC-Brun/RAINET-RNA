
library( ade4)
library(data.table)
# library(devtools)
# install_github("ggbiplot", "vqv")
require(ggplot2)
require(grid)
require(gridExtra)
library(ggbiplot)
source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

# testing
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/test/pca/interaction_score_matrix_altered2.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/mondal2015_ji2016/mondal2015_s5_ji2016/interaction_score_matrix.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/ribosomal/RNA_transport/fake_matrix_for_testing.tsv"

# all proteins
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/mondal2015_ji2016/annotated_interactions.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/cellular/sample_annotated_interactions.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/ribosomal/sample_annotated_interactions.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/biotypes/annotated_interactions.tsv"

# RBP only
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/cellular/RBP_only/annotated_interactions.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/ribosomal/RBP_only/annotated_interactions.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/biotypes/RBP_only/annotated_interactions.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/lncrna_vs_mrna/annotated_interactions.tsv"

#small datasets
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/snoRNA/annotated_interactions.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/cellular/RNA_transport/cytoplasmic_vs_nuclear.txt"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/ribosomal/RNA_transport/free_polysomal_subset.txt"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/cellular/exportins_importins/cytoplasmic_vs_nuclear.txt"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/cellular/exportins_importins_expression_data/cytoplasmic_vs_nuclear.txt"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/lncRNA_analysis/LncRNA_analysis/results/PCA_analysis/carlevaro2016/cellular/exportins_importins_expression_data/annotated_interactions.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

####################
## script options
####################

annotCol = "annotation" # name of the column carrying the annotation information
minimum_category_size = 10 # minimum number of items in category for it to be plotted #0 if no filter

####################
## initialisation
####################

# exclude small categories from analysis
datasetLess = dataset
for (i in unique(dataset$annotation) ){
  c = count(dataset$annotation == i)$freq[[2]]
  if (c < minimum_category_size){    datasetLess = datasetLess[datasetLess$annotation != i]  }
}
dataset = datasetLess

# use first column as dataframe row names
row.names(dataset) = dataset$RNAs
dataset$RNAs = NULL

# get list of proteins
proteins = c(names(dataset))
proteins = proteins[1:length(proteins)-1]
proteins = gsub("|", '.', proteins, fixed = 1)

# # backup annotation data, remove annotation column from dataset
# annotationData = dataset[[annotCol]]
# dataset[[annotCol]] = NULL
nrow(dataset)
ncol(dataset)

# to allow column slice, need to convert data.table to data.frame
dataset = data.frame(dataset)

####################
## Run PCA
####################

## ADE4
acp <- dudi.pca( dataset[, proteins], scannf = FALSE, nf = ncol(dataset), scale= FALSE)

## Produce colors for dataset
sample_list = unique( dataset$annotation)
simple_color_palette = rainbow( length( sample_list))
simple_colors = unlist( lapply( sample_list, function( x){
  return( simple_color_palette[ which( sample_list == x)])
}))

# Plot the content of information by axes
barplot( 100*acp$eig/sum( acp$eig),
         names.arg = seq(1,length(acp$eig),1),
         xlab="Axes",
         ylab="Percentage of information",
         cex.main=0.7)

# Plot biplot
par( mfrow= c(2,2))
scatter(acp, clab.row = 0, clab.col = 0, posieig = "none", clabel=0, xax=1, yax=2)
s.class(dfxy = acp$li, fac = factor( dataset$annotation), col = simple_colors, add.plot = TRUE, cstar = 0, xax=1, yax=2)
scatter(acp, clab.row = 0, clab.col = 0, posieig = "none", clabel=0, xax=1, yax=3)
s.class(dfxy = acp$li, fac = factor( dataset$annotation), col = simple_colors, add.plot = TRUE, cstar = 0, xax=1, yax=3)
scatter(acp, clab.row = 0, clab.col = 0, posieig = "none", clabel=0, xax=2, yax=3)
s.class(dfxy = acp$li, fac = factor( dataset$annotation), col = simple_colors, add.plot = TRUE, cstar = 0, xax=2, yax=3)
scatter(acp, clab.row = 0, clab.col = 0, posieig = "none", clabel=0, xax=3, yax=4)
s.class(dfxy = acp$li, fac = factor( dataset$annotation), col = simple_colors, add.plot = TRUE, cstar = 0, xax=3, yax=4)

# Compute the best variables on the first 3 axes (code from Lionel Spinelli)

threshold_norm = 0.50
threshold_percent = 0.01
names_best_on_axis_all_axes = vector()
for( cs_index in c(1,2,3)){
  # Compute on the axis, the best variables, i.e. the variable which have the max coordinates on that axis
  current_vector = acp$c1[[cs_index]]
  number_best = floor( length( current_vector) * threshold_percent)
  decreasing_order = order( abs(current_vector), decreasing = TRUE)
  names_best_on_axis = row.names( acp$c1)[decreasing_order][1: number_best]
  cat("\n\nOn CS", cs_index,",", threshold_percent*100, "% of best variables are\n")
  print( names_best_on_axis)
  names_best_on_axis_all_axes = unique( append( names_best_on_axis_all_axes, names_best_on_axis))
  
  # Look at the variables with the maximal cumulative information on each axes
  cat("\nOn CS", cs_index, ", direct and cummulative information carried by best variables:\n")
  old_norm = 0
  for( gene_index in seq(1, length( current_vector))){
    current_norm = sqrt (sum( (acp$c1[[ cs_index]][decreasing_order][1: gene_index])^2))
    cat("\n", row.names( acp$c1)[decreasing_order][gene_index], "\t:\t", signif((current_norm - old_norm)*100, 3), "\t:\t", signif( current_norm*100 , 3), " %")
    if( current_norm >= threshold_norm){
      break
    }
    old_norm = current_norm
  }
}

## R prcomp
pca_result <- prcomp(dataset[, proteins], center = TRUE, scale. = FALSE)
# scale false since all data have same units/scale, they all come from catRAPID

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

listPlot = list()
listPlot[[1]] = ggbiplot(pca_result, choices = 1:2, var.axes = 0, var.scale = 1, alpha = 0.2, ellipse = TRUE,  groups = dataset$annotation) +
 scale_color_discrete(name = '') +
 theme(legend.direction = 'horizontal', legend.position = 'top')
listPlot[[2]] = ggbiplot(pca_result, choices = c(1,3), var.axes = 0, var.scale = 1, alpha = 0.2, ellipse = TRUE,  groups = dataset$annotation) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
listPlot[[3]] = ggbiplot(pca_result, choices = 2:3, var.axes = 0, var.scale = 1, alpha = 0.2, ellipse = TRUE,  groups = dataset$annotation) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
listPlot[[4]] = ggbiplot(pca_result, choices = c(3,4), var.axes = 0, var.scale = 1, alpha = 0.2, ellipse = TRUE,  groups = dataset$annotation) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

grid_arrange_shared_legend(listPlot)


# analysing the principal components

rot = pca_result$rotation
length(rot[,1])
tail(sort(rot[,1]))
tail(sort(rot[,2]))

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
