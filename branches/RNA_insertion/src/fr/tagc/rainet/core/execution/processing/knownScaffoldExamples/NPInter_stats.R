
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(Epi)
library(pROC)
library(ROCR)


args <- commandArgs(TRUE)
# Test if we have enough arguments
if( length(args) != 1){
  stop("Rscript: Bad argument number")
}
inputFile = args[1]

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/NPInterPredictionValidation/scores.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/NPInterPredictionValidation/scores_all_interactions.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Get only the validated for drawing points on them
validated = dataset[dataset$in_validated_set == 1,]
nonValidated = dataset[dataset$in_validated_set == 0]

print (ks.test(validated$catrapid_score, nonValidated$catrapid_score, alternative = c("two.sided")) )

plt0 <- ggplot(dataset, aes(x=catrapid_score, colour = as.factor(in_validated_set))) + 
  geom_density() + 
  theme_minimal() +
  xlab( "Score") +
  ylab( "Density") +
  annotate("text",  x=Inf, y = Inf, label = paste("# True: ",nrow(validated),"\n# False: ",nrow(dataset)-nrow(validated)), vjust=1, hjust=1)
plt0



### ROC Curve ###

calcThreshold <- function(response, predict) {
  r <- pROC::roc(response, predict)
  r$thresholds[which.max(r$sensitivities + r$specificities)]
}

# # Sub-sampling the non-validated dataset

repetitions = 100
nsubsample = nrow(validated)
print (nsubsample)

colours = sample(colours(), repetitions)

cutoffs = c()
aurocs = c()
kss = c()

for (i in seq(1, repetitions)){
  nonValidatedSample = nonValidated[sample( nrow(nonValidated), nsubsample ), ]
  datasetSample = rbind(nonValidatedSample, validated) 
  
  #   # Kolmogorow-Smirnov  
  #   ks = ks.test(validated$catrapid_score, nonValidatedSample$catrapid_score, alternative = c("two.sided"))$p.value
  #   kss = c(kss, ks)
  
  # Optimal cutoff / cutpoint
  cutoff = calcThreshold(datasetSample$in_validated_set, datasetSample$catrapid_score)
  cutoffs = c(cutoffs, cutoff)
  
  # ROC
  pred = prediction(datasetSample$catrapid_score, datasetSample$in_validated_set, label.ordering = NULL)
  perf = performance(pred, measure = "tpr", x.measure = "fpr")
  if (i == 1){
    plot(perf, lty=5, col= colours[i],main="")
  }else{
    plot(perf, lty=5, col= colours[i],main="", add=TRUE)
  }
  
  # AUROC
  auc = performance(pred, measure = "auc")
  aurocs = c(aurocs, auc@y.values[[1]])
  
}

# Means and stds of randomisations
mean(cutoffs)
sd(cutoffs)

mean(aurocs)
sd(aurocs)


# # Plotting ROC for the whole dataset (Time consuming)
# plt1 = ROC( form = in_validated_set ~ catrapid_score, data = dataset, plot="ROC", main="ROC curve", MI=TRUE, MX=TRUE, PV=TRUE) 
# print (calcThreshold(dataset$in_validated_set, dataset$catrapid_score))

png(filename=paste(inputFile,".png",sep=""))
dev.off()


### Getting max sensitivity OR max specificity

r <- pROC::roc(datasetSample$in_validated_set, datasetSample$catrapid_score)
#r$thresholds[which.max(r$sensitivities + r$specificities)]

#Note: -inf will always give max sensitivities
# one option is to set a "cost" using ROCR, other is to look for local maximas
r$thresholds[which.max(r$sensitivities)]
r$thresholds[which.max(r$specificities)]

# # local maximas

r$thresholds
print (coords(roc=r, x = "local maximas", ret='threshold') )

# # cost option
pred = prediction(datasetSample$catrapid_score, datasetSample$in_validated_set, label.ordering = NULL)

# increasing cost of false negatives, i.e. get cutoff where no false negative is allowed
cost = performance(pred, "cost", cost.fp = 1, cost.fn = 10000000)
pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
pred@fn[[1]][which.min(cost@y.values[[1]])]
pred@fp[[1]][which.min(cost@y.values[[1]])]

# increasing cost of false positives, i.e. get cutoff where no false positive is allowed
cost = performance(pred, "cost", cost.fp = 10000000, cost.fn = 1)
pred@cutoffs[[1]][which.min(cost@y.values[[1]])]
pred@fn[[1]][which.min(cost@y.values[[1]])]
pred@fp[[1]][which.min(cost@y.values[[1]])]


