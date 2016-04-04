
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(Epi)
library(pROC)

## This is script produces distribution of scores and ROC curve ##

args <- commandArgs(TRUE)
# Test if we have enough arguments
if( length(args) != 1){
  stop("Rscript: Bad argument number")
}
inputFile = args[1]

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/scoreColumn1/ENST00000483525_SAMMSON.tsv_scores.tsv"

dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

### Distribution ###

# Get only the validated for drawing points on them
validated = dataset[dataset$in_validated_set == "True",]


plt0 <- ggplot(dataset, aes(x=catrapid_score, colour = in_validated_set)) + 
  geom_density() + 
  theme_minimal() +
  xlab( "Score") +
  ylab( "Density") +
  annotate("text",  x=Inf, y = Inf, label = paste("# True: ",nrow(validated),"\n# False: ",nrow(dataset)-nrow(validated)), vjust=1, hjust=1)

plt1 <- ggplot(dataset, aes(x=catrapid_score)) + 
  geom_histogram(binwidth=0.05, colour="black", fill="white") +
  theme_minimal() +
  xlab( "Score") +
  ylab( "Protein count")
for (i in validated$catrapid_score){
  loop_input = paste("geom_point(aes(x=",i,",y=1, size=0.1, colour='red', position='dodge') )")
  plt1 <- plt1 + eval(parse(text=loop_input))  
}
plt1 <- plt1 + guides(colour=FALSE)

png(filename=paste(inputFile,".png",sep=""))

grid.arrange(plt0, plt1)
dev.off()


### ROC ###

dataset$in_validated_set = as.integer(as.logical(dataset$in_validated_set))

png(filename=paste(inputFile,"_ROC.png",sep=""))

plt2 = ROC( form = in_validated_set ~ catrapid_score, data = dataset, plot="ROC", main="ROC curve", MI=TRUE, MX=FALSE, PV=FALSE) 

calcThreshold <- function(response, predict) {
  r <- pROC::roc(response, predict)
  r$thresholds[which.max(r$sensitivities + r$specificities)]
}
# annotate with the catRAPID cutoff point / threshold
text(0.8,0.4, paste("Opt cutoff:", calcThreshold(dataset$in_validated_set, dataset$catrapid_score), sep = " " ) )

dev.off()
