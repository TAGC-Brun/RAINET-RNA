
library(ggplot2)
library(reshape)
require(grid)
require(gridExtra)
library(data.table)

# tx_expression_avg_file1 = "/home/diogo/testing/tx_expression_avg_muscle.tsv"
# tx_expression_avg_file2 = "/home/diogo/testing/tx_expression_avg_muscle_ipr.tsv"
# tx_expression_avg_file1 = "/home/diogo/testing/ProcessGTExData_testing/transcript_expression_metrics.tsv"
# tx_expression_avg_file2 = "/home/diogo/testing/ProcessGTExData_testing/transcript_expression_metrics_no_outliers.tsv"

# use fread for faster file reading than read.table
expression_df1 <- fread(tx_expression_avg_file, stringsAsFactors = FALSE, header = TRUE, sep="\t")
expression_df2 <- fread(tx_expression_avg_no_outliers_file, stringsAsFactors = FALSE, header = TRUE, sep="\t")

## Coefficient variation of RPKM

# Histogram of coefficients of variation before outlier removal
plt1 <- ggplot(expression_df1[expression_df1$coef_variation != 0,], aes(x = coef_variation ) )  +
  geom_histogram(binwidth = 0.1, fill="white", colour="black") + 
  labs(title="Distribution of coefficient variations across transcripts, including outliers") +
  labs(x="Coefficient variation", y="Frequency") + 
  theme_minimal()

# Histogram of coefficients of variation after outlier removal
plt2 <- ggplot(expression_df2[expression_df2$coef_variation != 0,], aes(x = coef_variation ) )  +
  geom_histogram(binwidth = 0.1, fill="white", colour="black") + 
  labs(title="Distribution of coefficient variations across transcripts, excluding outliers") +
  labs(x="Coefficient variation", y="Frequency") + 
  theme_minimal()

grid.arrange( plt1, plt2)


## RPKM means

print ("#### Summary without removing outliers ####")
print("summary(expression_df1$rpkm_mean)")
print(summary(expression_df1$rpkm_mean))

print ("#### Summary removing outliers ####")
print("summary(expression_df2$rpkm_mean)")
print(summary(expression_df2$rpkm_mean))


XLIM_MAX = 5
print(paste("Plots with xlim:",XLIM_MAX, sep=" "))

# Histogrames of mean RPKMs of all transcripts, before outlier removal
plt3 <- ggplot(expression_df1, aes(x = rpkm_mean ) )  +
  geom_histogram(binwidth = 0.1, fill="white", colour="black") + 
  geom_vline(aes(xintercept = mean( rpkm_mean)),color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept = 0.1),color="red", linetype="dashed", size=1) +
  labs(title="Mean transcript expression, including outliers") +
  labs(x="RPKM mean", y="Count") + 
  xlim(c(-0.1,XLIM_MAX)) + 
  theme_minimal()

# Histogrames of mean RPKMs of all transcripts, after outlier removal
plt4 <- ggplot(expression_df2, aes(x = rpkm_mean ) )  +
  geom_histogram(binwidth = 0.1, fill="white", colour="black") + 
  geom_vline(aes(xintercept = mean( rpkm_mean)),color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept = 0.1),color="red", linetype="dashed", size=1) +
  labs(title="Mean transcript expression, excluding outliers") +
  labs(x="RPKM mean", y="Count") + 
  xlim(c(-0.1,XLIM_MAX)) + 
  theme_minimal()

grid.arrange( plt3, plt4)


# plt2 <- ggplot(expression_df, aes(x = ExprMedian ) )  +
#   geom_histogram(binwidth = 0.1) + 
#   geom_vline(aes(xintercept=mean(ExprMedian)),color="blue", linetype="dashed", size=1) +
#   geom_vline(aes(xintercept=0.1),color="red", linetype="dashed", size=1) +
#   labs(title="Median transcript expression") +
#   labs(x="RPKM median", y="Count") + 
#   xlim(c(-0.1,XLIM_MAX)) + 
#   theme_minimal()
# plt2
# 
# # Density plots for all RPKMs of each tissue
# plt3 <- ggplot(expression_df, aes(x = rpkm_mean ) )  +
#   geom_line(stat="density",aes(color=TissueName)) +
#   geom_vline(aes(xintercept=0.1),color="red", linetype="dashed", size=1) +
#   labs(title="Mean transcript expression") +
#   labs(x="RPKM mean") + 
#   xlim(c(0.1,XLIM_MAX)) + 
#   theme_minimal()
# plt3
# Scatterplot of Mean vs coef_variation
#plt6 <- ggplot(expression_df, aes(x=rpkm_mean, y=coef_variation)) +
#  geom_point(shape=1,alpha=1/4)
#  geom_smooth(method=lm)
#  xlim(c(0.0,300)) 
#plt6



