
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

#### Report on data presence ####

inputFile = paste( output_folder, report_rna_expression_data_presence, sep = "/")

#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/real_linc_expr/expr1.0/Report/rna_expression.tsv"

inputData <- read.table(inputFile, header = TRUE, sep = "\t")

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.7)),
  colhead = list(fg_params=list(cex = 0.8)),
  rowhead = list(fg_params=list(cex = 0.8)))

myt <- gridExtra::tableGrob( inputData, theme = mytheme)

grid.draw(myt)


#### Report on data values ####

inputFile = paste( output_folder, report_rna_expression, sep = "/")

print (paste("INPUT FILE:",inputFile) )

dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

dataLimit = 0.8

print (paste("DISPLAYING DATA WITH XLIM AT PERCENTILE:",dataLimit) )

plt1 <- ggplot(dataset, aes(x = type, y = meanExpression, fill = type ) )  +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() + 
  scale_y_continuous(limits = quantile(dataset$meanExpression, c(0.0, dataLimit))) +
  theme_minimal()
plt1

plt2 <- ggplot(dataset, aes(x = transcriptBiotype, y = meanExpression, fill = type ) )  +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() + 
  scale_y_continuous(limits = quantile(dataset$meanExpression, c(0.0, dataLimit))) +
  theme_minimal()
plt2

grid.arrange( plt1, plt2)

# plt1 <- ggplot(dataset, aes(x = meanExpression ) )  +
#   geom_histogram( binwidth = 0.01) +
#   scale_x_continuous(limits =  c(0.0, 1.0)) +
#   ylab("# RNAs") +
#   xlab("Mean expression (RPKM)") +
#   theme_minimal()
# plt1


#### Report on proteins/rnas with expression per tissue ####

inputFile1 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/real_lncrna_expr/stats_full/Report/proteins_expressed_per_tissue.tsv"
inputFile2 = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/real_lncrna_expr/stats_full/Report/rnas_expressed_per_tissue.tsv"

dataset1 = read.table( inputFile1, stringsAsFactors = FALSE, header = TRUE, sep="\t")
dataset2 = read.table( inputFile2, stringsAsFactors = FALSE, header = TRUE, sep="\t")

plt1 <- ggplot(dataset1, aes(x = tissue, y = number_expressed, fill = "red" ) )  +
  geom_bar( ) +
  coord_flip() + 
  scale_fill_manual(values=c("#d7191c")) +
  ggtitle("Protein expression per tissue") +
  theme_minimal()
plt1

plt2 <- ggplot(dataset2, aes(x = tissue, y = number_expressed, fill = "blue" ) )  +
  geom_bar( ) +
  coord_flip() + 
  scale_fill_manual(values=c("#2b83ba")) +
  ggtitle("LncRNA expression per tissue") +
  theme_minimal()
plt2

dataset1$type = "Proteins"
dataset2$type = "RNAs"


# maxDataset1 = 15974 #max(dataset1$number_expressed)
# maxDataset2 = 500 #max(dataset2$number_expressed) # 10392 is the number of lincRNAs in GENCODE
# 
# dataset1$percentage_expressed = dataset1$number_expressed * 100 / maxDataset1
# dataset2$percentage_expressed = dataset2$number_expressed * 100 / maxDataset2
# 
# mergedDataset = rbind(dataset1,dataset2)
# 
# plt3 <- ggplot(mergedDataset, aes(x = tissue, y = percentage_expressed, fill = type ) )  +
#   geom_bar(position = "dodge", stat = "identity") +
#   coord_flip() + 
#   theme_minimal()
# plt3

# ###########################################################################################################################
# ## Plot for JOBIM 2016
# ###########################################################################################################################
# plt3 <- ggplot(mergedDataset, aes(x = tissue, y = percentage_expressed, fill = type ) )  +
#   geom_bar(position = "dodge", stat = "identity") +
#   coord_flip() + 
#   ylab("% molecules expressed") +
#   scale_fill_manual(values=c("#d7191c", "#2b83ba")) +
#   theme_minimal() + 
#   theme(text = element_text(size=11))
# plt3
# 
# 
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/analysisStrategy/real_linc_expr/expr1.0_rerun/Report/freq_number_tissue_expressed.txt"
# 
# dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
# 
# dataset$percentage = dataset$frequency * 100.0 / sum(dataset$frequency)
# 
# plt2 <- ggplot( dataset, aes(x = number_tissues, y = frequency)) + 
#   geom_bar( stat = "identity") + 
#   theme_minimal() +
#   theme(text = element_text(size=11)) + 
#   xlab( "Number of tissues with co-presence" ) +
#   ylab( "# interactions")
# plt2
# 
# sum( dataset$frequency[ dataset$number_tissues != 0])
# sum( dataset$percentage[ dataset$number_tissues != 0])
# 
