#!/usr/bin/env Rscript
# Useful R functions to be reused

require(grid)
require(gridExtra)

give.n <- function(x){   return(c(y = -10, label = length(x))) }

# Function to run tests for all pairwise comparisons
# arg dataset : a dataframe containing at least two columns, a categorical and a numeric
# arg metricToUse : column name of numerical data
# arg annotationCol : column name of categorical data
all_vs_all_tests <- function( dataset, metricToUse, annotationCol) { 
  
  # options for outputing values in scientific notation
  options( scipen = 0 )
  options( digits = 2 )

  # get all different categories  
  categories = unique(dataset[[annotationCol]])
  
  print (categories)

  # loop pairs of categories
  comparisons = c()
  doneComparisons = c()
  pvalues = c()
  means1 = c()
  means2 = c()
  
  for (i in categories){
    for (j in categories) {

        invComparison = paste(j, " vs ", i) %in% comparisons

        if (i != j && !invComparison ){
          comparison = paste(i, " vs ", j)
          print (comparison)
          set1 = dataset[ dataset[[annotationCol]] == i][[metricToUse]]
          set2 = dataset[ dataset[[annotationCol]] == j][[metricToUse]]
          means1 = c( means1, mean(set1))
          means2 = c( means2, mean(set2))
          pvalue = ks.test(set1, set2, alternative = c("two.sided"))$p.value
          comparisons = c(comparisons, comparison)
          pvalues = c( pvalues, pvalue)
          
          #keeping track of already performed comparisons
          doneComparisons = c(doneComparisons, comparison)
      }
    }
  }
  
  # correct pvalues
  correctedPvalues = p.adjust(pvalues, method = "bonferroni", n = length(pvalues))

  # create output table  
  statsDF = data.frame( comparisons, means1, means2, pvalues, correctedPvalues)
  grid.table( statsDF)
}
