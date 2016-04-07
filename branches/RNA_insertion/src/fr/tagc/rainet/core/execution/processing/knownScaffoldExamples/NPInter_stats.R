
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(Epi)
#library(pROC)
library(ROCR)

args <- commandArgs(TRUE)
# Test if we have enough arguments
if( length(args) != 1){
  stop("Rscript: Bad argument number")
}
inputFile = args[1]


#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/NPInterPredictionValidation/test_scores.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/NPInterPredictionValidation/scores.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/NPInterPredictionValidation/scores_all_interactions.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Separate validated from non-validated 
validated = dataset[dataset$in_validated_set == 1,]
nonValidated = dataset[dataset$in_validated_set == 0]

###################################### 
###################################### 
### Plot catRAPID scores distributions in and outside NPInter 
###################################### 
###################################### 

# testing difference of distributions (Kolmogorow-Smirnov  )
print (ks.test(validated$catrapid_score, nonValidated$catrapid_score, alternative = c("two.sided")) )

plt0 <- ggplot(dataset, aes(x=catrapid_score, colour = as.factor(in_validated_set))) + 
  geom_density() + 
  theme_minimal() +
  xlab( "Score") +
  ylab( "Density") +
  annotate("text",  x=Inf, y = Inf, label = paste("# True: ",nrow(validated),"\n# False: ",nrow(dataset)-nrow(validated)), vjust=1, hjust=1)
plt0


###################################### 
###################################### 
### ROC / cutoff analysis
###################################### 
###################################### 

###################################### 
###################################### 
## Random subsamples 
###################################### 
###################################### 

# Approach: Sub-sampling the non-validated dataset to match length of validated dataset, and keeping the validated dataset as it is

repetitions = 100
nsubsample = nrow(validated)
print (nsubsample)

colours = sample(colours(), repetitions)

cutoffs = c()
aurocs = c()
kss = c()

# Create ROC plots and other metrics
for (i in seq(1, repetitions)){
  nonValidatedSample = nonValidated[sample( nrow(nonValidated), nsubsample ), ]
  datasetSample = rbind(nonValidatedSample, validated) 
  
  #   # Kolmogorow-Smirnov  
  #   ks = ks.test(validated$catrapid_score, nonValidatedSample$catrapid_score, alternative = c("two.sided"))$p.value
  #   kss = c(kss, ks)
  
  # ROC using ROCR
  pred = prediction(datasetSample$catrapid_score, datasetSample$in_validated_set, label.ordering = NULL)
  perf = performance(pred, measure = "tpr", x.measure = "fpr")
  if (i == 1){
    plot(perf, lty=5, col= colours[i],main="", ylab = "Sensitivity", xlab = "1-Specificity")
  }else{
    plot(perf, lty=5, col= colours[i],main="", add=TRUE)
  }

  # Optimal cutoff / cutpoint
  cutoff = perf@alpha.values[[1]][which.max(1-perf@x.values[[1]] + perf@y.values[[1]])] 
  cutoffs = c(cutoffs, cutoff)
  
  # AUROC
  auc = performance(pred, measure = "auc")
  aurocs = c(aurocs, auc@y.values[[1]])
  
}

# Mean cutoff and std of randomisations
paste("Mean cutoff:",round(mean(cutoffs),2),"std:",round(sd(cutoffs),2))
# Mean auc and std of randomisations
paste("Mean auroc:",round(mean(aurocs),2),"std:",round(sd(aurocs),2))

###################################### 
# FDR curves
###################################### 
# Note: I'm in fact resampling, not using the same samples in previous plot.
# Problem is that I cannot store "plot" to a variable, and so, inside a loop I can only add data to one plot
# I tried storing FDR prediction results in a vector and then loop them to create new plot but does not work
# I cannot use performance objects to ggplots, I would have to create new data structure for ggplots
# I could potentially store the sampling data and loop them
for (i in seq(1, repetitions)){
  nonValidatedSample = nonValidated[sample( nrow(nonValidated), nsubsample ), ]
  datasetSample = rbind(nonValidatedSample, validated) 
  
  pred = prediction(datasetSample$catrapid_score, datasetSample$in_validated_set, label.ordering = NULL)
  perf = performance(pred, measure = "pcfall")
  if (i == 1){
    plot(perf, lty=5, col= colours[i],main="",ylab = "False discovery rate FP/(TP+FP)", xlim = c(min(dataset$catrapid_score),max(dataset$catrapid_score)), ylim = c(0,1))
  }else{
    plot(perf, lty=5, col= colours[i],main="", add=TRUE)
  }
}
abline(v = mean(cutoffs))


###################################### 
###################################### 
#### Whole dataset 
###################################### 
###################################### 

pred = prediction(dataset$catrapid_score, dataset$in_validated_set, label.ordering = NULL)

# ROC curve with ROCR package (fast)
plot(performance(pred, measure = "tpr", x.measure = "fpr"))

# Sensitivity plot
plot(performance(pred, measure = "sens"))

# Specificity plot
plot(performance(pred, measure = "spec"))

# False discovery rate on whole dataset
#pcfall: Prediction-conditioned fallout. P(Y = - | Yhat = +). Estimated as: FP/(TP+FP).
plot(performance(pred, measure = "pcfall"), ylab = "False discovery rate FP/(TP+FP)", ylim = c(0,1))

###################################### 
# Adding noise to dataset (missing information / stability analysis) 
###################################### 
# How does the ROC behave if I introduce some previously negatives as being positives

# Separate validated from non-validated 
validated = dataset[dataset$in_validated_set == 1,]
nonValidated = dataset[dataset$in_validated_set == 0]

repetitions = 10
sampleSize = nrow(validated) # how many negatives becoming positives

cutoffs = c()
aurocs = c()

for (i in seq(1, repetitions)){
  
  print(paste("Repetition:",i) )
  
  # randomly pick indexes from negatives of where to make change
  indexesToChange = sample( nrow(nonValidated), sampleSize )

  # create new copy of nonValidated, which will be modified
  copyNonValidated = nonValidated
  
  # now the non validated set will contain some validated items
  copyNonValidated$in_validated_set[indexesToChange] = 1

  # join 'non'-validated and validated
  datasetSample = rbind(validated,copyNonValidated) 
  
  # ROC using ROCR
  pred = prediction(datasetSample$catrapid_score, datasetSample$in_validated_set, label.ordering = NULL)
  perf = performance(pred, measure = "tpr", x.measure = "fpr")
  if (i == 1){
    plot(perf, lty=5, col= colours[i],main="", ylab = "Sensitivity", xlab = "1-Specificity")
  }else{
    plot(perf, lty=5, col= colours[i],main="", add=TRUE)
  }
  
  # Optimal cutoff / cutpoint
  cutoff = perf@alpha.values[[1]][which.max(1-perf@x.values[[1]] + perf@y.values[[1]])] 
  cutoffs = c(cutoffs, cutoff)
  
  # AUROC
  auc = performance(pred, measure = "auc")
  aurocs = c(aurocs, auc@y.values[[1]])
  
}
# Add the original (before noise introduced) line
pred = prediction(dataset$catrapid_score, dataset$in_validated_set, label.ordering = NULL)
perf = performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, lwd = 3, lty=1, col= "black",main="", add=TRUE)


# Mean cutoff and std of randomisations
paste("Mean cutoff:",round(mean(cutoffs),2),"std:",round(sd(cutoffs),2))
# Mean auc and std of randomisations
paste("Mean auroc:",round(mean(aurocs),2),"std:",round(sd(aurocs),2))




# png(filename=paste(inputFile,".png",sep=""))
# dev.off()


# ### Getting max sensitivity OR max specificity
# 
# r <- pROC::roc(datasetSample$in_validated_set, datasetSample$catrapid_score)
# #r$thresholds[which.max(r$sensitivities + r$specificities)]
# 
# 
# #Note: -inf will always give max sensitivities
# # one option is to set a "cost" using ROCR, other is to look for local maximas
# r$thresholds[which.max(r$sensitivities)]
# r$thresholds[which.max(r$specificities)]
# 
# # # local maximas
# 
# r$thresholds
# print (coords(roc=r, x = "local maximas", ret='threshold') )
# 
# # # cost option
# pred = prediction(datasetSample$catrapid_score, datasetSample$in_validated_set, label.ordering = NULL)
# 
# # increasing cost of false negatives, i.e. get cutoff where no false negative is allowed
# cost = performance(pred, "cost", cost.fp = 1, cost.fn = 10000000)
# pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
# pred@fn[[1]][which.min(cost@y.values[[1]])]
# pred@fp[[1]][which.min(cost@y.values[[1]])]
# 
# # increasing cost of false positives, i.e. get cutoff where no false positive is allowed
# cost = performance(pred, "cost", cost.fp = 10000000, cost.fn = 1)
# pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
# pred@fn[[1]][which.min(cost@y.values[[1]])]
# pred@fp[[1]][which.min(cost@y.values[[1]])]


# # Manually calculate FDR
# cutoff = -50
# 
# above = datasetSample[datasetSample$catrapid_score >= cutoff]
# below = datasetSample[datasetSample$catrapid_score < cutoff]
# nrow(above)
# nrow(below)
# 
# validated = datasetSample[datasetSample$in_validated_set == 1,]
# nonValidated = datasetSample[datasetSample$in_validated_set == 0]
# nrow(validated)
# nrow(nonValidated)
# 
# aboveValidated = above[above$in_validated_set == 1,]
# nrow(aboveValidated)
# 
# TP = nrow(aboveValidated) #nrow(validated)
# FP = nrow(above) - nrow(aboveValidated)
# 
# TP
# FP
# 
# FP / (FP + TP) 

##
# calcThreshold(datasetSample$in_validated_set, datasetSample$catrapid_score)
# calcThreshold <- function(response, predict) {
#   r <- pROC::roc(response, predict)
#   r$thresholds[which.max(r$sensitivities + r$specificities)]
# }

# # Plotting ROC for the whole dataset with Epi package (Time consuming)
# plt1 = ROC( form = in_validated_set ~ catrapid_score, data = dataset, plot="ROC", main="ROC curve", MI=TRUE, MX=TRUE, PV=TRUE) 
# # Calculate cutoff maximizing sensitivity and specificity
# print (calcThreshold(dataset$in_validated_set, dataset$catrapid_score))
