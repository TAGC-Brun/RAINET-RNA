
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)

inputFile = "/home/diogo/Documents/RAINET_data/catRAPID/catRAPID_all_vs_all/mrna/3utr/3utr_vs_mrna_correlation/score_correlation_shuf1000000.txt"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

median( dataset$ScoreFile1)
median( dataset$ScoreFile2)

correlation = cor(dataset$ScoreFile1, dataset$ScoreFile2, method = "spearman", use="complete")
correlationSign = as.numeric(cor.test(dataset$ScoreFile1, dataset$ScoreFile2, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")

plt1 <- ggplot(data = dataset, aes(x = ScoreFile1, y = ScoreFile2) )  +
  geom_point( shape = 1, alpha=1/4 ) +
  geom_smooth( ) +
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
  theme_minimal()
plt1


### for Gian ###

riccardoLibrary = c(97.84,35.8,118.99,80.55,45.15,57.35,114.34,49.36,-51.63,-17.81,4.07,26.57,78.38,-2.51,60.87,33.69,92.98,12.05,1.27,13.67)
diogoLibrary = c(62.46,27.6,136.69,107.27,41.11,60.53,83.35,29.7,-31.11,-10.94,3.58,24.58,61.07,10.02,49.26,30.14,75.01,6.52,-13.54,-1.15)
webserverLibrary = c(62.666,27.798,119.516,93.365,39.547,59.255,78.205,25.524,-29.81,-8.965,3.163,24.008,72.52,8.469,51.675,31.588,71.762,5.058,-11.618,-2.143)

dataset = data.frame(riccardoLibrary,diogoLibrary,webserverLibrary)


## Riccardo vs Diogo

correlation = cor(dataset$riccardoLibrary, dataset$diogoLibrary, method = "spearman", use="complete")
correlationSign = as.numeric(cor.test(dataset$riccardoLibrary, dataset$diogoLibrary, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")
plt1 <- ggplot(data = dataset, aes(x = riccardoLibrary, y = diogoLibrary) )  +
  geom_point( shape = 1, alpha=1/4 ) +
  geom_smooth( ) +
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
  theme_minimal()
plt1

## Diogo vs Webserver

correlation = cor(dataset$webserverLibrary, dataset$diogoLibrary, method = "spearman", use="complete")
correlationSign = as.numeric(cor.test(dataset$webserverLibrary, dataset$diogoLibrary, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")
plt2 <- ggplot(data = dataset, aes(x = webserverLibrary, y = diogoLibrary) )  +
  geom_point( shape = 1, alpha=1/4 ) +
  geom_smooth( ) +
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
  theme_minimal()
plt2

## Riccardo vs Webserver

correlation = cor(dataset$webserverLibrary, dataset$riccardoLibrary, method = "spearman", use="complete")
correlationSign = as.numeric(cor.test(dataset$webserverLibrary, dataset$riccardoLibrary, method = "spearman")$p.value)
correlationText = paste("Corr:", round(correlation,2),"(pval:", round(correlationSign),")")
plt3 <- ggplot(data = dataset, aes(x = webserverLibrary, y = riccardoLibrary) )  +
  geom_point( shape = 1, alpha=1/4 ) +
  geom_smooth( ) +
  annotate("text", x = Inf, y = Inf, label = correlationText, hjust = 1, vjust =1  ) +
  theme_minimal()
plt3


grid.arrange(plt1,plt2,plt3)
