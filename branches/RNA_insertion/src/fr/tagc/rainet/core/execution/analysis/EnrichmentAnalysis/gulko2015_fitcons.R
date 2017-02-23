
# 23-Feb-2017 Diogo Ribeiro
# Script to plot Coverage by annotation type in regions that exceed a fitCons score threshold of S (Gulko2015).

library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
#source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

lincRNALists = "/home/diogo/Documents/RAINET_data/lncRNA_info/gulko2013_fitCons/test/Rtest/lincRNA_plus_CDS_annotation.txt"
gulkoOutput = "/home/diogo/Documents/RAINET_data/lncRNA_info/gulko2013_fitCons/test/Rtest/tx_outfile_linc_CDS.txt"

lincRNAAnnotation <- fread(lincRNALists, stringsAsFactors = FALSE, header = TRUE, sep="\t")
gulkoData <- fread(gulkoOutput, stringsAsFactors = FALSE, header = TRUE, sep="\t")

mergeDataset <- merge( gulkoData, lincRNAAnnotation, by=c("transcript"))

## If wanting to add a "other" category 
#mergeDataset <- merge( gulkoData, lincRNAAnnotation, by=c("transcript"), all = TRUE)
#mergeDataset$group[is.na(mergeDataset$group)] = "other"

# for each group of transcripts and for each score threshold, calculate proportion of nucleotides above cutoff
counter = 0
for (grp in unique(mergeDataset$group) ){
  dtGroup = mergeDataset[ mergeDataset$group == grp]
  for (scr in unique(mergeDataset$score_cutoff)){
    counter = counter + 1
    dtGroupScore = dtGroup[dtGroup$score_cutoff == scr]
    percAbove = round(sum( dtGroupScore$ncAbove) / sum( dtGroupScore$txSize), 2)
    datalist[[counter]] <- c(grp,scr,percAbove)
   }
}

results = do.call(rbind.data.frame, datalist)
colnames(results) <- c("group", "score_threshold", "proportion_above")


## Bar plot with proportion per threshold
plt1 <- ggplot( results, aes(x = score_threshold, y = proportion_above, fill = group) )  +
  geom_bar( stat = "identity", position = "dodge") +
  theme_minimal()
plt1

results


## Bar plot with proportion per threshold
plt1 <- ggplot( results, aes(x = score_threshold, y = proportion_above, color = group) )  +
  geom_point( ) +
#  geom_line() +
  theme_minimal()
plt1

