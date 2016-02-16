
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(gridExtra)


#From: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

expression_input_file = "/home/diogo/testing/expression_test.csv"
#expression_input_file = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/dataFiles/expression.csv"

nc <- max(count.fields(expression_input_file, sep=","))
expression_df <- read.table(expression_input_file, sep=",", row.names = 1, col.names=paste("V",1:nc,sep="."), fill=T)
expression_df <- as.data.frame(t(expression_df))
expression_df <- melt(expression_df)
colourCount = length(unique(expression_df$variable))
getPalette = colorRampPalette(brewer.pal(7, "Set1")) #this creates a function, which will take total number of shades of colours to use

tissues = factor(expression_df$variable)

plt1 <- ggplot(expression_df, aes(x = value ) )  +
  geom_density(aes(color=tissues), size=3) +
  theme_minimal()

plt2 <- ggplot(expression_df, aes(x = value ) )  +
  geom_density(aes(fill=factor(variable)), size=3) +
  theme_minimal()

grid_arrange_shared_legend(plt1,plt2)


