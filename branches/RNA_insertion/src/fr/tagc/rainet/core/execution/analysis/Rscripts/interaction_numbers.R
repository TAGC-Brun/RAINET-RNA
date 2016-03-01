
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/real/Report/interaction_numbers.tsv"
#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_expected/Report/interaction_numbers.tsv"
#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_expected/Report_w_filter/interaction_numbers.tsv"

inputFile = paste( output_folder, report_interaction_numbers, sep = "/")

print (paste("INPUT FILE:",inputFile) )

inputData <- read.table(inputFile, header = TRUE, sep = "\t")

# The values for total interactions and RNA numbers are in very different scales, I will plot them separately

# #
# Total interactions
totalInteractions <- inputData[,1:2]
totalInteractions <- melt(totalInteractions)

# Bar plot for numbers of items before and after filter
plt1 <- ggplot( totalInteractions, aes(x = variable, y = value, fill = Data)) + 
  geom_bar( stat = "identity", position=position_dodge()) + #aes(colour = Data)
  coord_flip() + 
  theme_minimal() +
  geom_text(aes(label = value), hjust = 0.5, size = 4, position = position_dodge( width=1)) +
  scale_fill_discrete( breaks = factor(totalInteractions$Data, levels = rev(levels(totalInteractions$Data)))) +
  xlab( "Number of interactions") +
  ylab( "Frequency")

# #
# RNA numbers
rnaNumbers <- inputData[,c(1,3:ncol(inputData) )] #header plus all values except total interactions
rnaNumbers <- melt(rnaNumbers)

# Re-ordering labels
rnaNumbers$variable <- factor(rnaNumbers$variable, levels = rev(levels(rnaNumbers$variable)))

# Bar plot for numbers of items before and after filter
plt2 <- ggplot( rnaNumbers, aes(x = variable, y = value, fill = Data)) + 
  geom_bar( stat = "identity", position=position_dodge()) + #aes(colour = Data)
  coord_flip() + 
  theme_minimal() +
  geom_text(aes(label = value), hjust = 0.5, size = 4, position = position_dodge( width=1)) +
  scale_fill_discrete( breaks = factor(rnaNumbers$Data, levels = rev(levels(rnaNumbers$Data)))) +
  xlab( "Biotype") +
  ylab( "Frequency")


# Percentages plot

proportion <- inputData[2,2:ncol(inputData)]*100.0 / inputData[1,2:ncol(inputData)]
proportion <- melt(proportion)

proportion$variable <- factor(proportion$variable, levels = rev(levels(proportion$variable)))

# Bar plot for numbers of items before and after filter
plt3 <- ggplot( proportion, aes(x = variable, y = value)) + 
  geom_bar( stat = "identity", position=position_dodge(), fill = "red") +
  coord_flip() + 
  geom_text( aes( label = round(value,2) ) ) + #, hjust = -0.5, size = 4, position = position_dodge( width=1), color = "white") +
  theme_minimal() +
  ylim(c(0.0,100.0)) + 
  xlab( "Biotype") +
  ylab( "Percentage after interaction filter")

grid.arrange( plt1, plt2, plt3)
