
# 23-Feb-2017 Diogo Ribeiro
# Script to plot Coverage by annotation type in regions that exceed a fitCons score threshold of S (Gulko2015).

library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
#source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

rnaLists = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/Gulko2015/transcript_annotation/transcript_groups.txt"
# gulkoOutput = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/Gulko2015/read_gulko.out"
gulkoOutput = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/Gulko2015/read_gulko_large_range.out"

rnaAnnotation <- fread(rnaLists, stringsAsFactors = FALSE, header = TRUE, sep="\t")
gulkoData <- fread(gulkoOutput, stringsAsFactors = FALSE, header = TRUE, sep="\t")

mergeDataset <- merge( gulkoData, rnaAnnotation, by=c("transcript"))

# nrow(rnaAnnotation[rnaAnnotation$group == "CDS"]) == nrow(mergeDataset[mergeDataset$group == "CDS"])/14
nrow(rnaAnnotation[rnaAnnotation$group == "CDS"]) == nrow(mergeDataset[mergeDataset$group == "CDS"])/76
# 76 is the number of different thresholds used

## If wanting to add a "other" category 
#mergeDataset <- merge( gulkoData, rnaAnnotation, by=c("transcript"), all = TRUE)
#mergeDataset$group[is.na(mergeDataset$group)] = "other"

# for each group of transcripts and for each score threshold, calculate proportion of nucleotides above cutoff
counter = 0
datalist = list()
for (grp in unique(mergeDataset$group) ){
  dtGroup = mergeDataset[ mergeDataset$group == grp]
  for (scr in unique(mergeDataset$score_cutoff)){
    counter = counter + 1
    dtGroupScore = dtGroup[dtGroup$score_cutoff == scr]
    propAbove = round(sum( dtGroupScore$ncAbove) / sum( dtGroupScore$txSize), 2)
    datalist[[counter]] <- c(grp,scr,propAbove)
   }
}

results = do.call(rbind.data.frame, datalist)
colnames(results) <- c("group", "score_cutoff", "proportion_above")

# need to convert to numeric otherwise they are factor
results$score_cutoff = as.numeric(levels(results$score_cutoff)[results$score_cutoff]) 
results$proportion_above = as.numeric(levels(results$proportion_above)[results$proportion_above])

# ## Density plot with proportion per transcript
# plt0 <- ggplot( mergeDataset[mergeDataset$score_cutoff == 0.5], aes(x = percAbove, color = group) )  +
#   geom_density( ) +
#   theme_minimal()
# plt0

# ## Bar plot with proportion per threshold
# plt1 <- ggplot( results, aes(x = score_cutoff, y = proportion_above, fill = group) )  +
#   geom_bar( stat = "identity", position = "dodge") +
#   theme_minimal()
# plt1


## Scatter plot with proportion per threshold
plt2 <- ggplot( results, aes(x = score_cutoff, y = proportion_above, color = group) )  +
  geom_point( size = 3) +
  geom_line( size = 2) +
  xlab("fitCons score threshold") +
  ylab("Sequence covered") +
#  scale_x_continuous(breaks = c(0.06, 0.07, 0.08, 0.09, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80)) +
  scale_x_continuous(breaks = seq(from = 0, to = 0.8, by = 0.05))  +
  scale_colour_brewer(palette = "Set1", labels = c("3'UTR","5'UTR","CDS","Scaffolding lincRNAs","Other lincRNAs")) +
  labs(color = "") +
  #  scale_color_manual(values = simple_color_palette) +
  theme_bw() + 
  theme(legend.position="right", legend.direction= "vertical", legend.key = NULL) + 
  theme(text = element_text(size=25))
plt2

#print as 25.91 vs 5.98 inches

