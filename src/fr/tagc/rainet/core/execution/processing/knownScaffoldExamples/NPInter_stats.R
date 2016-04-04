
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)

args <- commandArgs(TRUE)
# Test if we have enough arguments
if( length(args) != 1){
  stop("Rscript: Bad argument number")
}
inputFile = args[1]
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/NPInterPredictionValidation/scores.tsv"

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

library(Epi)
library(pROC)

# datasetSample <- dataset[sample( nrow(dataset), 100000), ]

plt1 = ROC( form = in_validated_set ~ catrapid_score, data = dataset, plot="ROC", main="ROC curve", MI=TRUE, MX=TRUE, PV=FALSE) 

calcThreshold <- function(response, predict) {
  r <- pROC::roc(response, predict)
  r$thresholds[which.max(r$sensitivities + r$specificities)]
}
print (calcThreshold(dataset$in_validated_set, dataset$catrapid_score) )


png(filename=paste(inputFile,".png",sep=""))

dev.off()

