
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/real/Report/rna_numbers.tsv"
inputFile = paste( output_folder, report_rna_numbers, sep = "/")

print (paste("INPUT FILE:",inputFile) )

inputData <- read.table(inputFile, header = TRUE, sep = "\t")

rnaNumbers <- melt(inputData)

# Re-ordering labels
rnaNumbers$variable <- factor(rnaNumbers$variable, levels = rev(levels(rnaNumbers$variable)))

# Bar plot for numbers of items before and after filter
plt1 <- ggplot( rnaNumbers, aes(x = variable, y = value, fill = Data)) + 
  geom_bar( stat = "identity", position=position_dodge()) + #aes(colour = Data)
  coord_flip() + 
  theme_minimal() +
  geom_text(aes(label = value), hjust = 1, size = 4, position = position_dodge( width=1)) +
  scale_fill_discrete( breaks = factor(rnaNumbers$Data, levels = rev(levels(rnaNumbers$Data)))) +
  xlab( "Class of RNA") +
  ylab( "Frequency")

print(plt1)
