#!/usr/bin/env Rscript
# Useful R functions to be reused

require(grid)
require(gridExtra)

give.n <- function(x){   return(c(y = -10, label = length(x))) }

# Function to run Smirkov-Kolmorov tests for all pairwise comparisons
# arg dataset : a dataframe containing at least two columns, a categorical and a numeric
# arg metricToUse : column name of numerical data
# arg annotationCol : column name of categorical data
all_vs_all_tests <- function( dataset, metricToUse, annotationCol, verbose = 1) { 
  
  # options for outputing values in scientific notation
  options( scipen = 0 )
  options( digits = 2 )

  # get all different categories  
  categories = unique(dataset[[annotationCol]])

  # loop pairs of categories
  comparisons = c()
  pvalues = c()
  means1 = c()
  means2 = c()
  
  for (i in categories){
    for (j in categories) {

        invComparison = paste(j, " vs ", i) %in% comparisons
        
        if (i != j && !invComparison ){
          comparison = paste(i, " vs ", j)
          set1 = dataset[ dataset[[annotationCol]] == i][[metricToUse]]
          set2 = dataset[ dataset[[annotationCol]] == j][[metricToUse]]
          means1 = c( means1, mean(set1, na.rm=TRUE))
          means2 = c( means2, mean(set2, na.rm=TRUE))
          pvalue = ks.test(set1, set2, alternative = c("two.sided"))$p.value
          comparisons = c(comparisons, comparison)
          pvalues = c( pvalues, pvalue)
      }
    }
  }
  
  # correct pvalues
  correctedPvalues = p.adjust(pvalues, method = "bonferroni", n = length(pvalues))

  # create output table  
  statsDF = data.frame( comparisons, means1, means2, pvalues, correctedPvalues)
  if (verbose){
    grid.table( statsDF)
  }
  return (statsDF)
}




# Specific function to run plots for low complexity
# arg dataset1 : data frame, should have a column called annotation
# arg dataset2 : data frame, should have a column called annotation
# arg annotationCol : column name of categorical data
# arg metricToUse1 : column name of numerical data
plot_filt_dataset <- function( dataset1, dataset2, filtAnnot, metricToUse1) { 
  annotCol = "with_low_complexity"
  
  # filter dataset by type of annotation  
  if (filtAnnot == "*"){
    # no filter
    filtDataset1 = dataset1
    filtDataset2 = dataset2 }
  else{
    filtDataset1 = dataset1[dataset1$annotation == filtAnnot]
    filtDataset2 = dataset2[dataset2$annotation == filtAnnot]  }
  
#   # separate RBPs with or without low complexity region
  filtDataset1$with_low_complexity = filtDataset1$perc_low > 0
  filtDataset2$with_low_complexity = filtDataset2$perc_low > 0
#   filtDataset1$with_low_complexity = filtDataset1$total_length > 200 #500
#   filtDataset2$with_low_complexity = filtDataset2$total_length > 200 #500
  
  # table with number of items  
  table1 = count(filtDataset1, annotCol)
  names(table1) = c(annotCol,filtAnnot)

  # run external function for all vs all tests
  table2 = all_vs_all_tests( filtDataset1, metricToUse1, annotCol, verbose = 0)
  table3 = all_vs_all_tests( filtDataset2, metricToUse1, annotCol, verbose = 0)
  
  grid.arrange(
    tableGrob(table1),
    tableGrob(table2),
    tableGrob(table3),
    nrow=3)
  
  # compare distributions between values of with_low_complexity for the two given datasets 
  plt3 <- ggplot(data = filtDataset1, aes_string(x = metricToUse1, colour = "with_low_complexity") )  +
    geom_density( size = 1 ) +
    ggtitle( paste(filtAnnot, filtDataset1$type) ) +
    theme_minimal()
  
  plt4 <- ggplot(data = filtDataset2, aes_string(x = metricToUse1, colour = "with_low_complexity") )  +
    geom_density( size = 1 ) +
    ggtitle( paste(filtAnnot, filtDataset2$type) ) +
    theme_minimal()
  
  grid.arrange( plt3, plt4)
}



# ##########################################################################################
# Define a function that offers to automatically layout a list of plots in a single viewport
multiplot <- function( plots=NULL, ncols=1, nrows) {
  require(grid)
  
  numplots = length( plots)
  if( numplots > ncols * nrows){
    stop( "multiplot: To many plots for the numbers of rows and columns:", nrows, "*", ncols)
  }
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport( layout = grid.layout( nrows, ncols)))
  
  # Make each plot, in the correct location
  col=1
  row=1
  for (i in 1:numplots) {
    print(plots[[i]], vp = viewport(layout.pos.row = row,
                                    layout.pos.col = col))
    col = col + 1
    if( col > ncols){
      col = 1
      row = row + 1
    }
  }
}


#EDITED From template: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(plots) {
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}
