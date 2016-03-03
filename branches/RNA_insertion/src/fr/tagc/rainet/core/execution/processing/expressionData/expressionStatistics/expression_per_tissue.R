
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(gridExtra)

#expression_input_file = "/home/diogo/testing/expression_test.csv"
#expression_input_file = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/dataFiles/expression.csv"
#expression_input_file = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/dataFiles/expression_sample.csv"

nc <- max(count.fields(expression_input_file, sep=","))

expression_df <- read.table(expression_input_file, sep=",", row.names = 1, col.names=paste("V",1:nc,sep="."), fill=T)
expression_df <- as.data.frame(t(expression_df))
expression_df <- melt(expression_df)

colourCount = length(unique(expression_df$variable))
getPalette = colorRampPalette(brewer.pal(7, "Set1")) #this creates a function, which will take total number of shades of colours to use

print("ALL SAMPLES: summary(expression_df$value[expression_df$variable == Muscle − Skeletal]")
print(summary(expression_df$value[expression_df$variable == "Muscle − Skeletal"]) )
print("ALL SAMPLES: summary(expression_df$value[expression_df$variable == Fallopian Tube]")
print(summary(expression_df$value[expression_df$variable == "Fallopian Tube"]) )

plt1 <- ggplot(expression_df, aes(x = variable, y = log2(value) ) )  +
  geom_boxplot(fill=getPalette(colourCount)) +
  coord_flip() + 
#  theme(axis.text=element_text(size=4)) +
  theme_minimal() +
  theme(text = element_text(size=10)) +
  xlab("Tissues")
#    scale_y_discrete(limits=c("0"))

# ## Same but for the other file
# expression_df <- read.table(expression_sample_input_file, sep=",", row.names = 1, col.names=paste("V",1:nc,sep="."), fill=T)
# expression_df <- as.data.frame(t(expression_df))
# expression_df <- melt(expression_df)
# plt2 <- ggplot(expression_df, aes(x = variable, y = log2(value) ) )  +
#   geom_boxplot(fill=getPalette(colourCount)) +
#   coord_flip() + 
#   theme_minimal() +
#   theme(text = element_text(size=10)) +
#   xlab("Tissues")
# #    scale_y_discrete(limits=c("0"))
# 
# print("SAMPLED: summary(expression_df$value[expression_df$variable == Brain]")
# print(summary(expression_df$value[expression_df$variable == "Brain"])) 
# print("SAMPLED: summary(expression_df$value[expression_df$variable == Fallopian Tube]")
# print(summary(expression_df$value[expression_df$variable == "Fallopian Tube"]) )
# 
# 
# grid.arrange(plt1, plt2, ncol=2)

print( plt1)
