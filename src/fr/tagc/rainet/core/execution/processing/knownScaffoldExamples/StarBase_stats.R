
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)
library(Epi)
library(ROCR)

##### Based on NPInter_stats.R

args <- commandArgs(TRUE)
# Test if we have enough arguments
if( length(args) != 2){
  stop("Rscript: Bad argument number")
}
inputFile = args[1]
outFolder = args[2]

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/lncRNAs/scores.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/mRNAs/scores.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/mRNAs/stringent/scores.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/StarBasePredictionValidation/mRNAs/extraStringent/scores.tsv"

setwd(outFolder)

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Separate validated from non-validated 
validated = dataset[dataset$in_validated_set == 1,]
nonValidated = dataset[dataset$in_validated_set == 0]

###################################### 
###################################### 
#### Whole dataset 
###################################### 
###################################### 

###################################### 
# Plot catRAPID scores distributions in and outside StarBase
###################################### 

# testing difference of distributions (Kolmogorow-Smirnov  )
ksTestWhole = ks.test(validated$catrapid_score, nonValidated$catrapid_score, alternative = c("two.sided"))

plt0 <- ggplot(dataset, aes(x=catrapid_score, colour = as.factor(in_validated_set))) + 
  geom_density() + 
  theme_minimal() +
  xlab( "Score") +
  ylab( "Density") +
  annotate("text",  x=Inf, y = Inf, label = paste("# True: ",nrow(validated),"\n# False: ",nrow(dataset)-nrow(validated)), vjust=1, hjust=1) +
  annotate("text",  x=Inf, y = 0, label = paste("KS test p-val: ", round(ksTestWhole$p.value,2)), vjust=1, hjust=1)
plt0

###################################### 
# Plot catRAPID scores versus CLIP reads mapped in StarBase
###################################### 

correlation = cor(validated$clipReads,validated$catrapid_score, method = "spearman")
correlationSign = as.numeric(cor.test(validated$catrapid_score, validated$clipReads, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")

plt1 <- ggplot(validated, aes(x = clipReads, y = catrapid_score)) + 
  geom_point(shape=1) + 
  geom_smooth(method=lm) + 
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  )
plt1

###################################### 
# ROC curve 
###################################### 

# Use ROCR
pred = prediction(dataset$catrapid_score, dataset$in_validated_set, label.ordering = NULL)

perf = performance(pred, measure = "tpr", x.measure = "fpr")

# ROC curve with ROCR package (fast)
plot(perf, ylab = "Sensitivity", xlab = "1-Specificity", lwd= 3, main = "ROC curve of whole dataset")
abline(a = 0, b = 1, col = "gray60")

auc = performance(pred, measure = "auc")
text(0.9, 0.4, paste("AUC:",round(auc@y.values[[1]],2) ) )

# Optimal cutoff / cutpoint, best sensitivity and best specificity
cutoff = perf@alpha.values[[1]][which.max(1-perf@x.values[[1]] + perf@y.values[[1]])] 
cutoff
text(0.7,0.45, paste("Cutoff[Max(sens+spec)]:",round(cutoff,2) ) )

colours = rainbow(5)

# Cutoff cost analysis. Augmenting weight of sensitivity
print ("cost cutoff sens 1-spec")
for (cost in seq(1,5)){
  calc = which.max(1-perf@x.values[[1]] + cost*perf@y.values[[1]])
  cutoff = perf@alpha.values[[1]][calc]
  sens = perf@y.values[[1]][calc]
  spec = 1-perf@x.values[[1]][calc]
  print (paste(cost, cutoff , sens, spec) )
  points(x = 1-spec, y = sens, col = colours[cost], pch = 16)
}
legend("bottomright", paste("cost =", 1:5), col = colours[1:5], pch = 16, cex = 0.8, title = "Sensitivity")


###################################### 
###################################### 
#### Random subsamples 
###################################### 
###################################### 
# Approach: Sub-sampling the non-validated dataset to match length of validated dataset, and keeping the validated dataset as it is

# parameters
repetitions = 100
nsubsample = nrow(validated)

colours = sample(colours(), repetitions)

cutoffs = c()
aurocs = c()
kss = c()

# Create ROC plots and other metrics
for (i in seq(1, repetitions)){
  # create subsample of non-validated set
  nonValidatedSample = nonValidated[sample( nrow(nonValidated), nsubsample ), ]
  datasetSample = rbind(nonValidatedSample, validated) 
  
  # ROC using ROCR on subsample
  pred = prediction(datasetSample$catrapid_score, datasetSample$in_validated_set, label.ordering = NULL)
  perf = performance(pred, measure = "tpr", x.measure = "fpr")
  if (i == 1){
    plot(perf, lty=5, col= colours[i], ylab = "Sensitivity", xlab = "1-Specificity", lwd= 3, main = "ROC curve of randomised set")
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
cutoffText = paste("Mean cutoff:",round(mean(cutoffs),2),"std:",round(sd(cutoffs),2))
# Mean auc and std of randomisations
aucText = paste("Mean AUC:",round(mean(aurocs),2),"std:",round(sd(aurocs),2))

## adding more info to the plot
abline(a = 0, b = 1, col = "gray60")
text(0.7,0.3, paste("Subsample size:", nsubsample) )
text(0.7,0.2, paste("Number randomisations:", repetitions) )
text(0.7,0.1, cutoffText )
text(0.7,0.0, aucText )


###################################### 
# Adding noise to dataset (missing information / stability analysis) 
###################################### 
# How does the ROC behave if I introduce some previously negatives as being positives

# parameters
repetitions = 5
sampleSize = nrow(validated) #/10 # how many negatives becoming positives

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
    plot(perf, lty=5, col= colours[i], ylab = "Sensitivity", xlab = "1-Specificity", main = "ROC curve, noise addition")
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
# Add the original (before noise introduced) ROC line
pred = prediction(dataset$catrapid_score, dataset$in_validated_set, label.ordering = NULL)
perf = performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf, lwd = 3, lty=1, col= "black",main="", add=TRUE)

# Mean cutoff and std of randomisations
cutoffText = paste("Mean cutoff:",round(mean(cutoffs),2),"std:",round(sd(cutoffs),2))
# Mean auc and std of randomisations
aucText = paste("Mean AUC:",round(mean(aurocs),2),"std:",round(sd(aurocs),2))

## adding more info to the plot
abline(a = 0, b = 1, col = "gray60")
text(0.7,0.3, paste("# Noise items:", sampleSize) )
text(0.7,0.2, paste("Number randomisations:", repetitions) )
text(0.7,0.1, cutoffText )
text(0.7,0.0, aucText )

