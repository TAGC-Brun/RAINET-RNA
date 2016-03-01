
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

######## Interaction Partners #########

#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/real/Report/interaction_partners.tsv"

inputFile = paste( output_folder, report_interaction_partners, sep = "/")

print (paste("INPUT FILE:",inputFile) )

nc <- max(count.fields(inputFile, sep=","))

inputData <- read.table(inputFile, sep=",", row.names = 1, col.names=paste("V",1:nc,sep="."), fill=T)

interaction <- as.data.frame(t(inputData))
interaction <- melt(interaction)


# Density plot of interaction partners per biotype
plt1 <- ggplot( interaction, aes(x = value)) + 
  geom_histogram(binwidth = 1, fill="white", colour="black") + 
  theme_minimal() +
  xlab( "Interaction partners") +
  ylab( "Transcript count")

# Density plot of interaction partners per biotype
plt2 <- ggplot( interaction, aes(x = value)) + 
  geom_density( aes( color = variable), size = 1) +
  theme_minimal() +
  xlab( "Interaction partners")

grid.arrange( plt1, plt2)
