
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)

inputFile1 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/lncRNAs/scores.tsv"
inputFile2 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/mRNAs/scores.tsv"

dataset1 <- fread(inputFile1, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2 <- fread(inputFile2, stringsAsFactors = FALSE, header = TRUE, sep="\t")

nrow(dataset1)
nrow(dataset2)

###################################### 
# Sample sets for efficiency purposes
###################################### 

divSampled = 10

nsubsample1 = nrow(dataset1) / divSampled
nsubsample2 = nrow(dataset2) / divSampled

dataset1Sample = dataset1[sample( nrow(dataset1), nsubsample1 ), ]
dataset2Sample = dataset2[sample( nrow(dataset2), nsubsample2 ), ]

###################################### 
# Plot catRAPID scores distributions between mRNAs and lncRNAs
###################################### 

mean1 = mean(dataset1Sample$catrapid_score)
mean2 = mean(dataset2Sample$catrapid_score)
meanText = paste("Mean lncRNA score: ", round(mean1,2), "\nMean mRNA score:", round(mean2,2))

# testing difference of distributions (Kolmogorow-Smirnov  )
ksTest = ks.test(dataset1Sample$catrapid_score, dataset2Sample$catrapid_score, alternative = c("two.sided"))

plt0 <- ggplot() + 
  geom_density(data = dataset1Sample, aes(x=catrapid_score, colour = "lncRNAs" ) ) + 
  geom_density(data = dataset2Sample, aes(x=catrapid_score, colour = "mRNAs"  )) + 
  theme_minimal() +
  xlab( "Score") +
  ylab( "Density") +
  annotate("text",  x=Inf, y = 0, label = paste("KS test p-val: ", round(ksTest$p.value,2)), vjust=1, hjust=1) +
  annotate("text",  x=Inf, y = Inf, label = meanText, vjust=1, hjust=1)
plt0

