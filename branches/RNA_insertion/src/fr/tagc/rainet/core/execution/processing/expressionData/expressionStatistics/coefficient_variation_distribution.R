
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

tx_expression_avg_file1 = "/home/diogo/testing/tx_expression_avg_muscle.tsv"
tx_expression_avg_file2 = "/home/diogo/testing/tx_expression_avg_muscle_ipr.tsv"
#tx_expression_avg_file1 = "/home/diogo/testing/tx_expression_avg_all_tissues_all_tx.tsv"
#tx_expression_avg_file2 = "/home/diogo/testing/tx_expression_avg_all_tissues_all_tx_ipr.tsv"

expression_df1 = read.table( tx_expression_avg_file1, stringsAsFactors = FALSE, header = TRUE, sep="\t")
expression_df2 = read.table( tx_expression_avg_file2, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# Histogram of coefficients of variation
plt1 <- ggplot(expression_df1, aes(x = CoefVariation ) )  +
  geom_histogram(binwidth = 0.1, fill="white", colour="black") + 
  labs(title="Distribution of coefficient variations across transcripts, including outliers") +
  labs(x="Coefficient variation", y="Count") + 
  theme_minimal()


# Histogram of coefficients of variation
plt2 <- ggplot(expression_df2, aes(x = CoefVariation ) )  +
  geom_histogram(binwidth = 0.1, fill="white", colour="black") + 
  labs(title="Distribution of coefficient variations across transcripts, excluding outliers") +
  labs(x="Coefficient variation", y="Count") + 
  theme_minimal()


grid.arrange( plt1, plt2)


# Histogram of coefficients of variation
plt3 <- ggplot(expression_df1[expression_df1$CoefVariation != 0,], aes(x = CoefVariation ) )  +
  geom_histogram(binwidth = 0.01, fill="white", colour="black") + 
  labs(title="Distribution of coefficient variations across transcripts, including outliers") +
  labs(x="Coefficient variation", y="Count") + 
#  ylim(0,100) +
  theme_minimal()


# Histogram of coefficients of variation
plt4 <- ggplot(expression_df2[expression_df2$CoefVariation != 0,], aes(x = CoefVariation ) )  +
  geom_histogram(binwidth = 0.001, fill="white", colour="black") + 
  labs(title="Distribution of coefficient variations across transcripts, excluding outliers") +
  labs(x="Coefficient variation", y="Count") + 
#  ylim(0,100) +
  theme_minimal()

#expression_df2[which(expression_df2$CoefVariation != 0),]

grid.arrange( plt3, plt4)

