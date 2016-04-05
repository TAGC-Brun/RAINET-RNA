
library(ggplot2)
library(gplots)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)
library(data.table)

## Produce 2D ranking plot for selected transcripts

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/scoreColumn1/ENST00000483525_SAMMSON.tsv_scores.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/scoreColumn1/ENST00000565493_NORAD.tsv_scores.tsv"

dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

datasetMaxOrder <- dataset[order (dataset$catrapid_score, decreasing = TRUE),]
datasetMeanOrder <- dataset[order( dataset$catrapid_score_mean, decreasing = TRUE),] 


maxRanking = c() # max ranking
meanRanking = c() # mean ranking
rssValues = c() # will store RSS, using same ranking as for loop
trueXXCoordinates = c()
trueYYCoordinates = c()
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

pointsDataset = data.frame( maxRanking, meanRanking, datasetMaxOrder$uniprotac, rssValues)

plt1 = ggplot(pointsDataset, aes( x = maxRanking, y = meanRanking)) +
  geom_point(aes( color = rssValues)) +
  scale_colour_gradient(low = "red", high = "cyan") +
  geom_smooth(method=lm)
for (i in seq(1:length(trueXXCoordinates))){
  loop_input = paste("geom_point(shape = 21, colour = 'black', fill = 'white', size = 3, aes(x=",trueXXCoordinates[i],",y=",trueYYCoordinates[i],") )")
  plt1 <- plt1 + eval(parse(text=loop_input))  
}
plt1


cor(pointsDataset$maxRanking, pointsDataset$meanRanking)

