
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

#tx_expression_input_folder = "/home/diogo/testing/folder"
#tx_expression_input_folder = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/dataFiles/transcript_expression"

#EDITED From template: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(plots) {
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

listPlot = list()
i = 0
for (file in list.files(tx_expression_input_folder)){
  expression_input_file = paste(tx_expression_input_folder,file,sep="/")
  
  nc <- max(count.fields(expression_input_file, sep=","))
  expression_df <- read.table(expression_input_file, sep=",", row.names = 1, col.names=paste("V",1:nc,sep="."), fill=T)
  expression_df <- as.data.frame(t(expression_df))
  expression_df <- melt(expression_df)
  colourCount = length(unique(expression_df$variable))
  getPalette = colorRampPalette(brewer.pal(7, "Set1")) #this creates a function, which will take total number of shades of colours to use
  
  tissues = factor(expression_df$variable)

  i=i+1
  
  listPlot[[i]] <- ggplot(expression_df, aes(x = value ) )  +
    geom_density(aes(color=tissues), size=1) + # could also use fill instead of color
    ggtitle(file) +
    theme_minimal()
}

grid_arrange_shared_legend(listPlot)

