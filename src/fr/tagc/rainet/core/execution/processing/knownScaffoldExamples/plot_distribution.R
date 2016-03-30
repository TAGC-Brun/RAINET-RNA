
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)

args <- commandArgs(TRUE)
# Test if we have enough arguments
if( length(args) != 1){
  stop("Rscript: Bad argument number")
}
inputFile = args[1]

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation/ENST00000565493_NORAD.tsv_zcores.tsv"

dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

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

