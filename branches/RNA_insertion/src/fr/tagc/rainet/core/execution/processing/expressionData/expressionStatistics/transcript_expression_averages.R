

library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

#tx_expression_avg_file = "/home/diogo/testing/tx_expression_avg.tsv"
#tx_expression_avg_file = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/dataFiles/tx_expression_avg.tsv"

expression_df = read.table( tx_expression_avg_file, stringsAsFactors = FALSE, header = TRUE, sep="\t")


print("summary(expression_df$ExprMean)")
print(summary(expression_df$ExprMean))
print("summary(expression_df$ExprStd)")
print(summary(expression_df$ExprStd))
print("summary(expression_df$ExprMedian)")
print(summary(expression_df$ExprMedian))

#NOTE THAT I USE XLIM! #IDEALLY I SHOULD SUM UP VALUES ABOVE XLIM
XLIM_MAX = 5
print(paste("Plots with xlim:",XLIM_MAX, sep=" "))

# Histogrames for all RPKMs of all transcripts
plt1 <- ggplot(expression_df, aes(x = ExprMean ) )  +
  geom_histogram(binwidth = 0.1, fill="white", colour="black") + 
  geom_vline(aes(xintercept=mean(ExprMean)),color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=0.1),color="red", linetype="dashed", size=1) +
  labs(title="Mean transcript expression") +
  labs(x="RPKM mean", y="Count") + 
  xlim(c(-0.1,XLIM_MAX)) + 
  theme_minimal()
plt1

plt2 <- ggplot(expression_df, aes(x = ExprMedian ) )  +
  geom_histogram(binwidth = 0.1) + 
  geom_vline(aes(xintercept=mean(ExprMedian)),color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=0.1),color="red", linetype="dashed", size=1) +
  labs(title="Median transcript expression") +
  labs(x="RPKM median", y="Count") + 
  xlim(c(-0.1,XLIM_MAX)) + 
  theme_minimal()
plt2

# Density plots for all RPKMs of each tissue
plt3 <- ggplot(expression_df, aes(x = ExprMean ) )  +
  geom_line(stat="density",aes(color=TissueName)) +
  geom_vline(aes(xintercept=0.1),color="red", linetype="dashed", size=1) +
  labs(title="Mean transcript expression") +
  labs(x="RPKM mean") + 
  xlim(c(0.1,XLIM_MAX)) + 
  theme_minimal()
plt3

plt4 <- ggplot(expression_df, aes(x = ExprMedian ) )  +
  geom_line(stat="density",aes(color=TissueName)) +
  geom_vline(aes(xintercept=0.1),color="red", linetype="dashed", size=1) +
  labs(title="Median transcript expression") +
  labs(x="RPKM median") + 
  xlim(c(0.1,XLIM_MAX)) + 
  guides(colour=FALSE) +
  theme_minimal()
plt4

grid.arrange( plt1, plt2)
grid.arrange( plt3, plt4)


# Scatterplot of Mean vs CoefVariation
#plt6 <- ggplot(expression_df, aes(x=ExprMean, y=CoefVariation)) +
#  geom_point(shape=1,alpha=1/4)
#  geom_smooth(method=lm)
#  xlim(c(0.0,300)) 
#plt6



