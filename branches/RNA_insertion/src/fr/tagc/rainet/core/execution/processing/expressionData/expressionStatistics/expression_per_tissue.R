
library(ggplot2)
library(reshape)
library(RColorBrewer)

#expression_input_file = "/home/diogo/testing/expression_test.csv"
#expression_input_file = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/dataFiles/expression.csv"

nc <- max(count.fields(expression_input_file, sep=","))

expression_df <- read.table(expression_input_file, sep=",", row.names = 1, col.names=paste("V",1:nc,sep="."), fill=T)
expression_df <- as.data.frame(t(expression_df))

#boxplot(expression_df)
expression_df <- melt(expression_df)

colourCount = length(unique(expression_df$variable))
getPalette = colorRampPalette(brewer.pal(7, "Set1")) #this creates a function, which will take total number of shades of colours to use


plt <- ggplot(expression_df, aes(x = variable, y = log2(value) ) )  +
  geom_boxplot(fill=getPalette(colourCount)) +
  coord_flip() + 
  theme_minimal() +
  xlab("Tissues")
#    scale_y_discrete(limits=c("0"))
print(plt)

