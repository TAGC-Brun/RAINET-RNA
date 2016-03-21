
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

inputFile = paste( output_folder, report_rna_numbers, sep = "/")
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/real/Report/rna_numbers.tsv"

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

#### Pie plots ####

# Examples from http://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization

# use only values after filtering
afterFilt <- rnaNumbers[ rnaNumbers$Data == "After_RNA_filter",]

# pie plot for broad categories of RNAs
broadCat <- afterFilt[ afterFilt$variable == "MRNA" | afterFilt$variable == "OtherRNA" | afterFilt$variable == "LncRNA",]
plt2 <- ggplot( broadCat, aes(x = "", y = value, fill = variable))  + 
  geom_bar(width = 1, stat = "identity") + #aes(colour = Data) +
  coord_polar(theta = "y") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    axis.text.x=element_blank()  ) +
  geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), label = value), size=5)


# pie plot for categories of lncRNAs
lncRNABiotypes <- c("protein_coding","X3prime_overlapping_ncrna","antisense","bidirectional_promoter_lncrna","known_ncrna","lincRNA","macro_lncRNA","processed_transcript","sense_intronic","sense_overlapping","TEC")
specCat <- afterFilt[ afterFilt$variable %in% lncRNABiotypes,]

# Bar plot for numbers of items before and after filter
plt3 <- ggplot( specCat, aes(x = "", y = value, fill = variable))  + 
  geom_bar(width = 1, stat = "identity") + #aes(colour = Data) +
  coord_polar(theta = "y") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    axis.text.x=element_blank()  ) +
  geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), label = value), size=5)

grid.arrange( plt2, plt3)


