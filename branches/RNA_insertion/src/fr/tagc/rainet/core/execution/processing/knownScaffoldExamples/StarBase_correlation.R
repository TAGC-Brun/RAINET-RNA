
library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
library(dplyr)

#source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

starbaseInteractions = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/lncRNAs/new_dataset/stats/starbase_interactions_reformatted.tsv"
catrapidPredictions = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/starbase_interactions/starbase_interactions/storedInteractions_reformatted.tsv"

starbaseData <- fread(starbaseInteractions, stringsAsFactors = FALSE, header = TRUE, sep="\t")
catrapidData <- fread(catrapidPredictions, stringsAsFactors = FALSE, header = TRUE, sep="\t")

starbaseData
catrapidData


# merge, keeping only itnersction
mergedData = merge(x = starbaseData, y = catrapidData, by = c("transcript","protein"))

mergedData = mergedData[ mergedData$clipReads != 0]


plt1 = ggplot(mergedData, aes( x = clipReads, y = score)) +
  geom_point( size = 1.5) +
  geom_smooth()
plt1


negatives = head(mergedData[ order(mergedData$clipReads)], 100)
positives = tail(mergedData[ order(mergedData$clipReads)], 100)

mean(negatives$score)
mean(positives$score)

#results are the opposite than expected..


## biocomplex

negatives = mergedData[ mergedData$bioComplex <= 2]
positives = mergedData[ mergedData$bioComplex >= 3]

mean(negatives$score)
mean(positives$score)

ks.test(negatives$score, positives$score)

# biocomplex 3 has higher scores than biocomplex 1

############################
# Expression data
############################

starbaseInteractions = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/lncRNAs/new_dataset/stats/starbase_interactions_reformatted.tsv"
# catrapidPredictions = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/lncRNAs/expression_dataset/interactions_expression_reformatted.tsv"
catrapidPredictions = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/expression_data/starbase_min100/storedInteractions_reformatted.tsv"

starbaseData <- fread(starbaseInteractions, stringsAsFactors = FALSE, header = TRUE, sep="\t")
catrapidData <- fread(catrapidPredictions, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# merge, keeping only itnersction
mergedData = merge(x = starbaseData, y = catrapidData, by = c("transcript","protein"))

mergedData = mergedData[ mergedData$clipReads != 0]

plt1 = ggplot(mergedData, aes( x = clipReads, y = score)) +
  geom_point( size = 1.5) +
  geom_smooth()
plt1

negatives = head(mergedData[ order(mergedData$clipReads)], 100)
positives = tail(mergedData[ order(mergedData$clipReads)], 100)

mean(negatives$score)
mean(positives$score)

ks.test(negatives$score, positives$score)

