
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

######### Interaction scores ##########

#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/interaction_scores.tsv"
inputFileScores = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/real/Report/interaction_scores.tsv"

inputDataScores <- read.table(inputFileScores, row.names = 1, header = FALSE, sep = ",", fill = TRUE)

interactionScores <- as.data.frame(t(inputDataScores))
interactionScores <- melt(interactionScores)

# Density plot of interaction partners per biotype
plt1 <- ggplot( interactionScores, aes(x = value)) + 
  geom_histogram(binwidth = 1, fill="white", colour="black") + 
  theme_minimal() +
  xlab( "Interaction scores") +
  ylab( "Transcript count")

# Density plot of interaction scores per biotype
plt2 <- ggplot( interactionScores, aes(x = value)) + 
  geom_density( aes( color = variable), size = 1) +
#  geom_line( stat = "density", aes( color = variable) ) +
  theme_minimal() +
  xlab( "Interaction scores")

grid.arrange( plt1, plt2)


######## Interaction Partners #########

#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/interaction_partners.tsv"
inputFilePartners = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/real/Report/interaction_partners.tsv"

inputDataPartners <- read.table(inputFile, row.names = 1, header = FALSE, sep = ",", fill = TRUE)

interactionPartners <- as.data.frame(t(inputDataPartners))
interactionPartners <- melt(interactionPartners)

# Density plot of interaction partners per biotype
plt3 <- ggplot( interactionPartners, aes(x = value)) + 
  geom_histogram(binwidth = 1, fill="white", colour="black") + 
  theme_minimal() +
  xlab( "Interaction partners") +
  ylab( "Transcript count")

# Density plot of interaction partners per biotype
plt4 <- ggplot( interactionPartners, aes(x = value)) + 
  geom_density( aes( color = variable), size = 1) +
  theme_minimal() +
  xlab( "Interaction partners")

grid.arrange( plt3, plt4)
