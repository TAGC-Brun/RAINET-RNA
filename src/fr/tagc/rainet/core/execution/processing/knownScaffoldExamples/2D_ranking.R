
library(ggplot2)
library(gplots)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)
library(data.table)

## Produce 2D ranking plot for selected transcripts

# # SAMMSON
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/scoreColumn1/ENST00000483525_SAMMSON.tsv_scores.tsv"
#  
# # NORAD
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/scoreColumn1/ENST00000565493_NORAD.tsv_scores.tsv"
# 
# # XIST
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/scoreColumn1/ENST00000429829_XIST_001.tsv_scores.tsv"
# 
# XIST Mouse
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/scoreColumn1/all_interactions.85891.XIST_MOUSE_RBP_DBP_DISORDER.txt_scores.tsv"
# 
# # MEG3
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/scoreColumn1/ENST00000423456_MEG3_002.tsv_scores.tsv"
# 
# # NEAT1
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/scoreColumn1/ENST00000501122_NEAT1_001.tsv_scores.tsv"

# # NORAD vs CatRAPID
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/catRAPID_NORAD_rank.tsv"

dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

datasetMaxOrder <- dataset[order (dataset$catrapid_score, decreasing = TRUE),]
datasetMeanOrder <- dataset[order( dataset$catrapid_score_mean, decreasing = TRUE),] 

## Loop through one of the rankings, get coordinate of the other ranking, add other data
maxRanking = c() # max ranking
meanRanking = c() # mean ranking
rssValues = c() # will store RSS, using same ranking as for loop
trueXXCoordinates = c() # to store coordinate of validated proteins
trueYYCoordinates = c() # to store coordinate of validated proteins
# loop over one of the ordered sets, retrieve ranking position in the other set
for (maxIndex in seq( nrow( datasetMaxOrder)) ) {
  protID = datasetMaxOrder$uniprotac[ maxIndex]
  meanIndex = match( protID, datasetMeanOrder$uniprotac)
  maxRanking = c(maxRanking, maxIndex)
  meanRanking = c(meanRanking, meanIndex)
  # rss = sqrt( datasetMaxOrder$catrapid_score[maxIndex]^2 + datasetMaxOrder$catrapid_score_mean[maxIndex]^2)
  rss = sqrt( maxIndex^2 + meanIndex^2)
  rssValues = c( rssValues, rss)

  if (datasetMaxOrder$in_validated_set[maxIndex] == "True"){
      trueXXCoordinates = c(trueXXCoordinates, maxIndex)
      trueYYCoordinates = c(trueYYCoordinates, meanIndex)
  }
}

# create new dataframe which connects the two rankings
pointsDataset = data.frame( maxRanking, meanRanking, datasetMaxOrder$uniprotac, rssValues)

fileName = strsplit(inputFile,split="/")[[1]]
fileName = tail(fileName, n=1)

## 2D ranking plot

# correlation, and its significance
correlation = cor(pointsDataset$maxRanking, pointsDataset$meanRanking, method = "spearman")
correlationSign = as.numeric(cor.test(dataset$catrapid_score, dataset$catrapid_score_mean, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")

# Scatterplot of ranking and location of validated interactions
plt1 = ggplot(pointsDataset, aes( x = maxRanking, y = meanRanking)) +
  geom_point(aes( color = rssValues)) +
  scale_colour_gradient(low = "red", high = "cyan") +
  geom_smooth(method=lm) + 
  annotate("text", x = max(pointsDataset$maxRanking), y = max(pointsDataset$meanRanking), label = correlationText, hjust=1 ) +
  ggtitle( fileName)
for (i in seq(1:length(trueXXCoordinates))){
  loop_input = paste("geom_point(shape = 21, colour = 'black', fill = 'white', size = 3, aes(x=",trueXXCoordinates[i],",y=",trueYYCoordinates[i],") )")
  plt1 <- plt1 + eval(parse(text=loop_input))  
}
plt1

## Scatterplot with the raw points

# correlation, and its significance
correlation = cor(dataset$catrapid_score, dataset$catrapid_score_mean, method = "spearman")
correlationSign = as.numeric(cor.test(dataset$catrapid_score, dataset$catrapid_score_mean, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")

plt2 = ggplot(dataset, aes( x = catrapid_score, y = catrapid_score_mean)) +
  geom_point() +
  geom_smooth(method=lm) + 
  annotate("text", x = max(dataset$catrapid_score), y = max(dataset$catrapid_score_mean), label = correlationText, hjust=1 ) +
  labs(x="max score",y="mean score") +
  ggtitle( fileName)
plt2

# #grid.arrange(plt2,plt1)
# 
# # same plot but with different colours depending of validated or not
# plt2 = ggplot(dataset, aes( x = catrapid_score, y = catrapid_score_mean, fill = in_validated_set)) +
#   geom_point(aes(fill = in_validated_set, colour = in_validated_set)) +
#   geom_smooth(method=lm) + 
#   annotate("text", x = max(dataset$catrapid_score), y = max(dataset$catrapid_score_mean), label = correlationText, hjust=1 ) +
#   labs(x="max score",y="mean score") +
#   ggtitle( fileName)
# plt2


# # # same plot but with different colours depending of validated or not
# correlation = cor(dataset$catrapid_score, dataset$catrapid_score_mean, method = "spearman")
# correlationSign = as.numeric(cor.test(dataset$catrapid_score, dataset$catrapid_score_mean, method = "spearman")$p.value)
# correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")
# 
# plt2 = ggplot(dataset, aes( x = catrapid_score, y = catrapid_score_mean)) +
#   geom_point() +
#   geom_smooth(method=lm) + 
#   annotate("text", x = max(dataset$catrapid_score), y = max(dataset$catrapid_score_mean), label = correlationText, hjust=1 ) +
#   labs(x="catRAPID score",y="NORAD score") +
#   ggtitle( fileName)
# plt2
# 
