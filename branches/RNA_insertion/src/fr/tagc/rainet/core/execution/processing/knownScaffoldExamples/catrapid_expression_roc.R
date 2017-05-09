
## 02-May-2016 Diogo Ribeiro
## Script to create a ROC curve of a predictor based on both catRAPID and expression data
## Used generally for all types of experimental datasets, NPInter, StarBase, eCLIP etc.

library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(Epi)
library(ROCR)

###################################### 
# Read input data
###################################### 

# catRAPID file has scores and also whether interaction is experimental or not
catRAPIDfile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/NPInterPredictionValidation/allInteractions/scores_reformatted3.tsv"
# expression file only has expression data per interaction
expressionFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/expression_data/NPInter_interactions_ensembl68/storedInteractions_reformatted.tsv"

catRAPIDData = fread(catRAPIDfile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
expressionData = fread(expressionFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

###################################### 
# Prepare input data
###################################### 

# # catrapid scores only
# mergedData = catRAPIDData
# mergedData$metric = mergedData$catrapid_score # this should be rerun after any expression data transformation

# # expression data only
# mergedData = merge(x = catRAPIDData, y = expressionData, by = c("transcript","protein"))
# mergedData$metric = mergedData$expression_score

# merge catrapid and expression data, keeping only intersection
# mergedData = merge(x = catRAPIDData, y = expressionData, by = c("transcript","protein"))
# write.table(mergedData, file='/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/eCLIPPredictionValidation/ROC/only_112_rbps_15230_gencode_basic_lncRNAs/scores_plus_expression.tsv', quote=FALSE, sep='\t', col.names = T, row.names = F)

# # Doing a left merge, keeping all catRAPID interactions
mergedData = merge(x = catRAPIDData, y = expressionData, by = c("transcript","protein"), all.x = TRUE)
# Adding 0 to expression value when it would be NA
mergedData[is.na(mergedData)] = 0

# # Apply log transformation
# mergedData$expression_score = mergedData$expression_score + 1
# mergedData$expression_score = log10(mergedData$expression_score)

# Apply [0,1] normalisation
mergedData$catrapid_score = (mergedData$catrapid_score-min(mergedData$catrapid_score))/(max(mergedData$catrapid_score)-min(mergedData$catrapid_score))
mergedData$expression_score = (mergedData$expression_score-min(mergedData$expression_score))/(max(mergedData$expression_score)-min(mergedData$expression_score))

# summary(mergedData$expression_score)
# plt1 = ggplot(mergedData, aes( x = expression_score)) +
#   geom_histogram()
# plt1

# plt1 = ggplot(positives, aes( x = catrapid_score, y = expression_score)) +
#   geom_point( size = 1.5) +
#   geom_smooth()
# plt1
# cor(mergedData$catrapid_score, mergedData$expression_score)

# Separate validated from non-validated
validated = mergedData[mergedData$in_validated_set == 1,]
nonValidated = mergedData[mergedData$in_validated_set == 0]

### Subsample the negatives, to the same number as positives
nonValidated = nonValidated[sample(nrow(nonValidated), nrow(validated))]
mergedData = rbind(validated, nonValidated)

####################
### Basic statistics
####################
paste("Number of experimental interactions:", nrow(validated))
paste("Number of non-experimental interactions:",nrow(nonValidated))
# On catRAPID scores
paste("Median score of experimental interactions:",median(validated$catrapid_score))
paste("Median score of non-experimental interactions:",median(nonValidated$catrapid_score))
paste("Kolmogorov-Smirnov pvalue:", ks.test(validated$catrapid_score, nonValidated$catrapid_score, alternative = c("two.sided"))$p.value)
# On expression
paste("Median score of experimental interactions:",median(validated$expression_score))
paste("Median score of non-experimental interactions:",median(nonValidated$expression_score))
paste("Kolmogorov-Smirnov pvalue:", ks.test(validated$expression_score, nonValidated$expression_score, alternative = c("two.sided"))$p.value)


###################################### 
# ROC curve 
###################################### 

betas = c(0,0.001,0.01,0.1,1,10,100,1000)

for (beta in betas){

  # Define the metric to be used for ROC curve
  mergedData$metric = mergedData$catrapid_score + mergedData$expression_score * beta

  # Use ROCR
  pred = prediction(mergedData$metric, mergedData$in_validated_set, label.ordering = NULL)
  perf = performance(pred, measure = "tpr", x.measure = "fpr")
  
  auc = performance(pred, measure = "auc")
  
  print(paste("Beta:", beta, "AUC:", round(auc@y.values[[1]],2)))

}


#########################################
# Drawing ROC curve
#########################################

##################
# composite catrapid+expression model
##################
beta = 200
mergedData$metric = mergedData$catrapid_score + mergedData$expression_score * beta

# Use ROCR
pred = prediction(mergedData$metric, mergedData$in_validated_set, label.ordering = NULL)
perf = performance(pred, measure = "tpr", x.measure = "fpr")
auc = performance(pred, measure = "auc")

# ROC curve with ROCR package (fast)
plot(perf, ylab = "Sensitivity", xlab = "1-Specificity", lwd= 3, col = "dark green")
abline(a = 0, b = 1, lty = 5, col = "gray60")

text(0.5, 0.9, paste("AUC:",round(auc@y.values[[1]],2) ) )

##################
# catrapid only model
##################

beta = 0
mergedData$metric = mergedData$catrapid_score + mergedData$expression_score * beta

pred = prediction(mergedData$metric, mergedData$in_validated_set, label.ordering = NULL)
perf = performance(pred, measure = "tpr", x.measure = "fpr")
auc = performance(pred, measure = "auc")

lines(perf@x.values[[1]], perf@y.values[[1]], lwd = 3, lty = 2, col = "brown")

text(0.5, 0.65, paste("AUC:",round(auc@y.values[[1]],2) ) )

text(0.15,0.90, paste("# Positives:",nrow(validated) ) )
text(0.16,0.85, paste("# Negatives:",nrow(nonValidated) ) )

legend(0.4, 0.2, legend = c("catRAPID + expression","catRAPID only"), col = c("dark green","brown"), lty=1:2, cex=0.9)

### other tested files
# catRAPIDfile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/eCLIPPredictionValidation/ROC/only_112_rbps/score_reformatted.tsv"
# catRAPIDfile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/eCLIPPredictionValidation/ROC/only_112_rbps_15230_gencode_basic_lncRNAs/scores_reformatted3.tsv"
# catRAPIDfile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/eCLIPPredictionValidation/ROC/eclip_read_count/eclip_HepG2_positives_negatives.out"
# catRAPIDfile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/eCLIPPredictionValidation/ROC/eclip_read_count/eclip_K562_positives_negatives.out"
# catRAPIDfile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl68/interactions_with_npinter_positives/scores_reformatted3.tsv"
# catRAPIDfile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/lncRNAs/scores_reformatted3.tsv"
# expressionFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/expression_data/eCLIP_interactions/storedInteractions_reformatted.tsv"

