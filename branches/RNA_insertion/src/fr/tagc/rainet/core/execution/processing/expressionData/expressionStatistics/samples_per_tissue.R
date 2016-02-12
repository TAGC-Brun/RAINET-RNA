
library(RColorBrewer)
library(ggplot2)

# Load the data
annotation_df = read.table( annotation_input_file, stringsAsFactors = FALSE, header = FALSE, sep=",")

sorted_annotation_df <- annotation_df[order(-annotation_df$V2),] 

darkcols <- brewer.pal(length(sorted_annotation_df$V1), "Dark2")

barplot(sorted_annotation_df$V2, main="Samples per tissue", ylab="# Samples", names.arg=sorted_annotation_df$V1, las = 2,  col=darkcols) 

#plt <- ggplot(data=sorted_annotation_df, aes(x=sorted_annotation_df$V1, y=sorted_annotation_df$V2)) +
#  geom_bar(stat="identity") +
#  ylab("# genes") +
#  xlab("Species")
#print(plt)

print(paste("Total # of tissues : ",length(sorted_annotation_df$V2)) )
print(paste("Total # of samples : ",sum(sorted_annotation_df$V2)) )
print(paste("Minimum # samples per tissue: ",min(sorted_annotation_df$V2)) )
print(paste("Maximum # samples per tissue: ",max(sorted_annotation_df$V2)) )
