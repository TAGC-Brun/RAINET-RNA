
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)


###### For Barcelona meeting ######


### RNA numbers based on ncRNA fasta file

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/temp_results/files_for_barcelona_meeting/custom_rna_numbers_from_ncrna_fasta.txt"

inputData <- read.table(inputFile, header = TRUE, sep = "\t")

rnaNumbers <- melt(inputData)

# Re-ordering labels
rnaNumbers$variable <- factor(rnaNumbers$variable, levels = rev(levels(rnaNumbers$variable)))

lncRNABiotypes <- c("3prime_overlapping_ncrna","antisense","lincRNA","processed_transcript","sense_intronic","sense_overlapping","TEC")
afterFilt <- rnaNumbers[ rnaNumbers$Data == "After_RNA_filter",]

specCat <- afterFilt[ afterFilt$variable %in% lncRNABiotypes,]

sum(specCat$value)

plt3 <- ggplot( specCat, aes(x = variable, y = value, fill = variable))  + 
  geom_bar(aes(colour = variable), stat = "identity", position=position_dodge()) + 
  geom_text(aes(label=value, y = mean(value)), hjust=1, color="black", position = position_dodge(0.9), size=4)+
  coord_flip() +
  theme_minimal()
plt3

### RNA lengths based on ncRNA fasta file

inputFile = "/home/diogo/Documents/RAINET_data/Ensembl/13_03_2016/fasta_lengths"

inputData <- read.table(inputFile, header = FALSE)
inputData <- melt(inputData)

# To make a final bin which is values above xlim, alter list of values in order to have them at same bin
# be careful to then not use that vector for calculating mean etc, just for the plot
XLIM = 6000
vector <- c()
i = 1
for (val in inputData$value){
  if (val > XLIM){
    vector[i] <- XLIM  
  } else {
    vector[i] <- val
  }
  i=i+1
}

length(vector)  
length(inputData$value)

summary(vector)
summary(inputData$value)

vector <- melt(vector)

# Density plot of interaction partners per biotype
plt1 <- ggplot( vector, aes(x = value)) + 
  geom_histogram(binwidth = 100, fill="white", colour="black") + 
  geom_vline(aes(xintercept = mean( inputData$value)),color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept = 1200),color="red", linetype="dashed", size=1) +
  theme_minimal() +
  xlab( "Transcript length") +
  ylab( "Transcript count")
#  xlim(c(-0.1,XLIM))
plt1


# how many lncRNAs above catRAPID cutoff
length(inputData$value[inputData$value < 1200])

25046 / 30887


##### plot with interaction cutoff filtering #####


inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/temp_results/files_for_barcelona_meeting/custom_interaction_numbers.txt"

inputData <- read.table(inputFile, header = TRUE, sep = "\t")

rnaNumbers <- melt(inputData)

# Re-ordering labels
rnaNumbers$variable <- factor(rnaNumbers$variable, levels = rev(levels(rnaNumbers$variable)))
rnaNumbers$data <- factor(rnaNumbers$data, levels = rev(c("No filter", ">=0", ">=25", ">=50") ) )

levels(rnaNumbers$data)

plt4 <- ggplot( rnaNumbers, aes(x = variable, y = value, fill = data))  + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_text(aes(label=value, y = mean(value)), hjust=1, color="black", position = position_dodge(0.9), size=4)+
  coord_flip() +
  theme_minimal()
plt4


