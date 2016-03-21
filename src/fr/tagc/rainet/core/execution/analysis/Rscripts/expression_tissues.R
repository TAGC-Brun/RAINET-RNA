
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)
library(plyr)

inputFile = paste( output_folder, report_interactions_tissues_where_expressed, sep = "/")

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/real/Report/interactions_tissues_where_expressed.tsv"

print (paste("INPUT FILE:",inputFile) )

dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

if (nrow(dataset) != 0 ){
  #### Bar plot with tissue presence frequencies ####
  
  # read list_of_tissues into a single vector
  tissueVector = c()
  for (entry in dataset$list_of_tissues){
    tissues = strsplit(entry, ",")
    for (tiss in tissues){
      tissueVector = c(tissueVector,tiss)
    }  
  }
  # produce count frequencies from a single vector
  tissueFrequencies = count(tissueVector)
  
  # Bar plot for numbers of items before and after filter
  plt1 <- ggplot( tissueFrequencies, aes(x = x, y = freq)) + 
    geom_bar( stat = "identity", position=position_dodge()) + #aes(colour = Data)
    coord_flip() + 
    theme_minimal() +
    xlab( "Tissues") +
    ylab( "Frequency")
  
  print(plt1)
  
  #### Histogram of number of tissues with co-presence ####
  
  dataset <- melt(dataset)
  
  summary(dataset$value)
  
  plt2 <- ggplot( dataset, aes(x = value)) + 
    geom_histogram(binwidth = 1, fill="black", colour="black") + 
    theme_minimal() +
    xlab( "Number of tissues with co-presence" ) +
    ylab( "Interaction count")
  
  print(plt2)
} else{
  print ("No data for this section.")
}


