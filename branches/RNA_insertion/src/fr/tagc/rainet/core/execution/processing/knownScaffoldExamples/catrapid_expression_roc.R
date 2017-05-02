
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

catRAPIDfile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/eCLIPPredictionValidation/ROC/only_112_rbps/score_reformatted.tsv"
expressionFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/ReadCatrapid/Ensembl82/lncrna/expression_data/eCLIP_interactions/storedInteractions_reformatted.tsv"

catRAPIDData = fread(catRAPIDfile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
expressionData = fread(expressionFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# merge, keeping only intersection
mergedData = merge(x = catRAPIDData, y = expressionData, by = c("transcript","protein"))

# mergedData$catrapid_score = 0

# # # need to remove zero expression before doing log10
mergedData = mergedData[ mergedData$expression_score > 0]
mergedData$expression_score = log10(mergedData$expression_score)

# plt1 = ggplot(positives, aes( x = catrapid_score, y = expression_score)) +
#   geom_point( size = 1.5) +
#   geom_smooth()
# plt1

cor(mergedData$catrapid_score, mergedData$expression_score)

# Separate validated from non-validated
validated = mergedData[mergedData$in_validated_set == 1,]
nonValidated = mergedData[mergedData$in_validated_set == 0]

### Basic statistics
paste("Number of experimental interactions:", nrow(validated))
paste("Number of non-experimental interactions:",nrow(nonValidated))
paste("Median score of experimental interactions:",median(validated$catrapid_score))
paste("Median score of non-experimental interactions:",median(nonValidated$catrapid_score))
paste("Kolmogorov-Smirnov pvalue:", ks.test(validated$catrapid_score, nonValidated$catrapid_score, alternative = c("two.sided"))$p.value)

### Basic statistics
paste("Median score of experimental interactions:",median(validated$expression_score))
paste("Median score of non-experimental interactions:",median(nonValidated$expression_score))
paste("Kolmogorov-Smirnov pvalue:", ks.test(validated$expression_score, nonValidated$expression_score, alternative = c("two.sided"))$p.value)




###################################### 
# ROC curve 
###################################### 

betaMin = 0
betaMax = 20
betaStep = 5

for (beta in seq(betaMin, betaMax, by = betaStep)){

  # Define the metric to be used for ROC curve
  mergedData$metric = mergedData$catrapid_score + mergedData$expression_score * beta
  
  # Use ROCR
  pred = prediction(mergedData$metric, mergedData$in_validated_set, label.ordering = NULL)
  perf = performance(pred, measure = "tpr", x.measure = "fpr")
  
  auc = performance(pred, measure = "auc")
  
  print(paste("Beta:", beta))
  print(paste("AUC:", round(auc@y.values[[1]],2)))
  
}


beta = 20
mergedData$metric = mergedData$catrapid_score + mergedData$expression_score * beta

# Use ROCR
pred = prediction(mergedData$metric, mergedData$in_validated_set, label.ordering = NULL)
perf = performance(pred, measure = "tpr", x.measure = "fpr")

auc = performance(pred, measure = "auc")

# ROC curve with ROCR package (fast)
plot(perf, ylab = "Sensitivity", xlab = "1-Specificity", lwd= 3, main = "ROC curve")
abline(a = 0, b = 1, col = "gray60")
text(0.9, 0.4, paste("AUC:",round(auc@y.values[[1]],2) ) )

# Optimal cutoff / cutpoint, best sensitivity and best specificity
cutoff = perf@alpha.values[[1]][which.max(1-perf@x.values[[1]] + perf@y.values[[1]])] 
text(0.8,0.45, paste("Max(sens+spec):",round(cutoff,2) ) )

text(0.16,0.90, paste("# True:",nrow(validated) ) )
text(0.2,0.85, paste("# False:",nrow(nonValidated) ) )

colours = rainbow(3)

# Cutoff cost analysis. Augmenting weight of sensitivity
print ("cost cutoff sens 1-spec")
for (cost in seq(1,3)){
  calc = which.max(1-perf@x.values[[1]] + cost*perf@y.values[[1]])
  cutoff = perf@alpha.values[[1]][calc]
  sens = perf@y.values[[1]][calc]
  spec = 1-perf@x.values[[1]][calc]
  #   print (paste(cost, cutoff , sens, spec) )
  points(x = 1-spec, y = sens, col = colours[cost], pch = 16)
}
legend("bottomright", paste("cost =", 1:3), col = colours[1:3], pch = 16, cex = 0.8, title = "Sensitivity")
