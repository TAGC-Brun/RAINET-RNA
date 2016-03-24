
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/ENST00000565493_NORAD.tsv_zcores.tsv"

dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Get only the validated for drawing points on them
validated = dataset[dataset$in_validated_set == "True",]

plt0 <- ggplot(dataset, aes(x=catrapid_zscore, colour = in_validated_set)) + 
  geom_density() + 
  theme_minimal() +
  xlab( "Z-score") +
  ylab( "Density") +
  annotate("text",  x=Inf, y = Inf, label = paste("# True: ",nrow(dataset)-nrow(validated),"\n# False: ",nrow(validated)), vjust=1, hjust=1)
plt0

plt1 <- ggplot(dataset, aes(x=catrapid_zscore)) + 
  geom_histogram(binwidth=0.05, colour="black", fill="white") +
  theme_minimal() +
  xlab( "Z-score") +
  ylab( "Protein count")
for (i in validated$catrapid_zscore){
  loop_input = paste("geom_point(aes(x=",i,",y=1, size=0.1, colour='red', position='dodge') )")
  plt1 <- plt1 + eval(parse(text=loop_input))  
  }
plt1

grid.arrange(plt0, plt1)
