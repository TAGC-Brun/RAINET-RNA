
library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)

#source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

# experimentalInteractions = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/lncRNAs/new_dataset/stats/starbase_interactions_reformatted.tsv"
# experimentalInteractions = "/home/diogo/Documents/RAINET_data/ENCODE/eCLIP_from_Alex_Armaos/HepG2/eclip_interactions_readcount.out"
experimentalInteractions = "/home/diogo/Documents/RAINET_data/ENCODE/eCLIP_from_Alex_Armaos/K562/eclip_interactions_readcount.out"

# catrapidPredictions = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/starbase_interactions/starbase_interactions/storedInteractions_reformatted.tsv"
# catrapidPredictions = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl68/interactions_with_ensembl68_proteins_only/storedInteractions_reformatted.tsv"
catrapidPredictions = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl68/otherRNAs_only/storedInteractions_reformatted.tsv"
# catrapidPredictions = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/eCLIP_interactions/112_RBPs_15230_gencode_basic_lncRNAs/storedInteractions_reformatted.tsv"

experimentalData <- fread(experimentalInteractions, stringsAsFactors = FALSE, header = TRUE, sep="\t")
catrapidData <- fread(catrapidPredictions, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# merge, keeping only itnersction
mergedData = merge(x = experimentalData, y = catrapidData, by = c("transcript","protein"))
mergedData = mergedData[ mergedData$clipReads != 0]

nrow(experimentalData)
nrow(catrapidData)
nrow(mergedData)

plt1 = ggplot(mergedData, aes( x = clipReads, y = score)) +
  geom_point( size = 1.5) +
  # xlim( c(0,2200)) +
  geom_smooth()
plt1

negatives = head(mergedData[ order(mergedData$clipReads)], 100)
positives = tail(mergedData[ order(mergedData$clipReads)], 100)
mean(negatives$score)
mean(positives$score)
ks.test(negatives$score, positives$score)

# write.table(positives, file='/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/eCLIPPredictionValidation/ROC/eclip_read_count/eclip_K562_positives_ensembl68.tsv', quote=FALSE, sep='\t', col.names = T, row.names = F)
# write.table(negatives, file='/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/eCLIPPredictionValidation/ROC/eclip_read_count/eclip_K562_negatives_ensembl68.tsv', quote=FALSE, sep='\t', col.names = T, row.names = F)
# write.table(positives, file='/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/eCLIPPredictionValidation/ROC/eclip_read_count/eclip_HepG2_positives.tsv', quote=FALSE, sep='\t', col.names = T, row.names = F)
# write.table(negatives, file='/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/eCLIPPredictionValidation/ROC/eclip_read_count/eclip_HepG2_negatives.tsv', quote=FALSE, sep='\t', col.names = T, row.names = F)

## measure correlation between transcript length and score
correlation = cor(mergedData$clipReads, mergedData$score, method = "pearson", use="complete")
correlationSign = as.numeric(cor.test(mergedData$clipReads, mergedData$score, method = "pearson", use='complete')$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")
correlationText

# ## biocomplex
# 
# negatives = mergedData[ mergedData$bioComplex <= 2]
# positives = mergedData[ mergedData$bioComplex >= 3]
# 
# mean(negatives$score)
# mean(positives$score)
# 
# ks.test(negatives$score, positives$score)
# 
# # biocomplex 3 has higher scores than biocomplex 1

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

