
library(RColorBrewer)
library(ggplot2)

#annotation_input_file = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/dataFiles/annotation.csv"
  
# Load the data
annotation_df = read.table( annotation_input_file, stringsAsFactors = FALSE, header = FALSE, sep=",")

sorted_annotation_df <- annotation_df[order(annotation_df$V2),] 

colourCount = length(unique(sorted_annotation_df$V1))
getPalette = colorRampPalette(brewer.pal(7, "Set1")) #this creates a function, which will take total number of shades of colours to use

plt <- ggplot(data=sorted_annotation_df, aes(x=sorted_annotation_df$V1, y=sorted_annotation_df$V2)) +
  geom_bar(stat="identity", fill=getPalette(colourCount), color="white") + #if adding colour (fill) here, I will not have colour legend
  geom_text(aes(label=sorted_annotation_df$V2), color=getPalette(colourCount), hjust = -0.1, vjust=0.5, size=3.5) +
  scale_x_discrete(limits = sorted_annotation_df$V1) + # so that order of items respects the same order as in the input and not alphabectic
  theme_minimal() +
  coord_flip() +
  ylab("# samples") +
  xlab("Tissue")
print(plt)

print(paste("Total # of tissues : ",length(sorted_annotation_df$V2)) )
print(paste("Total # of samples : ",sum(sorted_annotation_df$V2)) )
#print(paste("Minimum # samples per tissue: ",min(sorted_annotation_df$V2)) )
#print(paste("Maximum # samples per tissue: ",max(sorted_annotation_df$V2)) )
