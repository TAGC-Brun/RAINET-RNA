
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

######### Interaction scores ##########

inputFile = paste( output_folder, report_interaction_scores, sep = "/")

#testing
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/temp_results/files_for_barcelona_meeting/data_without_interaction_filter/interaction_scores.tsv"
#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_expected/Report/interaction_scores.tsv"

print (paste("INPUT FILE:",inputFile) )

nc <- max(count.fields(inputFile, sep=","))

inputData <- read.csv(inputFile, sep=",", row.names = 1, col.names=paste("V",1:nc,sep="."), fill=T)

interaction <- as.data.frame(t(inputData))
interaction <- melt(interaction)

# Density plot of interaction scores per biotype
plt1 <- ggplot( interaction, aes(x = value)) + 
  geom_histogram(binwidth = 1, fill="white", colour="black") + 
  theme_minimal() +
  xlab( "Interaction scores") +
  ylab( "Transcript count")

# Density plot of interaction scores per biotype
plt2 <- ggplot( interaction, aes(x = value)) + 
  geom_density( aes( color = variable), size = 1) +
  #  geom_line( stat = "density", aes( color = variable) ) +
  theme_minimal() +
  xlab( "Interaction scores")

grid.arrange( plt1, plt2)
