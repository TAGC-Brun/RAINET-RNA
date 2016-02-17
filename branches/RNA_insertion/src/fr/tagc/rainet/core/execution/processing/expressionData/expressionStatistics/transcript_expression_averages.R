

library(ggplot2)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)

#tx_expression_avg_file = "/home/diogo/testing/tx_expression_avg.tsv"

#From template: https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
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

expression_df = read.table( tx_expression_avg_file, stringsAsFactors = FALSE, header = TRUE, sep="\t")

print("summary(expression_df$ExprMean)")
print(summary(expression_df$ExprMean))
print("summary(expression_df$ExprStd)")
print(summary(expression_df$ExprStd))
print("summary(expression_df$ExprMedian)")
print(summary(expression_df$ExprMedian))

#NOTE THAT I USE XLIM! #I SHOULD SUM UP VALUES ABOVE XLIM
XLIM_MAX = 5
print(paste("Plots with xlim:",XLIM_MAX, sep=" "))

plt1 <- ggplot(expression_df, aes(x = ExprMean ) )  +
  geom_histogram(binwidth = 0.1, fill="white", colour="black") + 
  geom_vline(aes(xintercept=mean(ExprMean)),color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=0.1),color="red", linetype="dashed", size=1) +
  labs(title="Mean transcript expression") +
  labs(x="RPKM mean", y="Count") + 
  xlim(c(-0.1,XLIM_MAX)) + 
  theme_minimal()
plt1

plt2 <- ggplot(expression_df, aes(x = ExprMedian ) )  +
  geom_histogram(binwidth = 0.1) + 
  geom_vline(aes(xintercept=mean(ExprMedian)),color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=0.1),color="red", linetype="dashed", size=1) +
  labs(title="Median transcript expression") +
  labs(x="RPKM median", y="Count") + 
  xlim(c(-0.1,XLIM_MAX)) + 
  theme_minimal()
plt2

plt3 <- ggplot(expression_df, aes(x = ExprMean ) )  +
  geom_line(stat="density",aes(color=TissueName)) +
  geom_vline(aes(xintercept=0.1),color="red", linetype="dashed", size=1) +
  labs(title="Mean transcript expression") +
  labs(x="RPKM mean") + 
  xlim(c(-0.1,XLIM_MAX)) + 
  theme_minimal()
plt3

plt4 <- ggplot(expression_df, aes(x = ExprMedian ) )  +
  geom_line(stat="density",aes(color=TissueName)) +
  geom_vline(aes(xintercept=0.1),color="red", linetype="dashed", size=1) +
  labs(title="Median transcript expression") +
  labs(x="RPKM median") + 
  xlim(c(-0.1,XLIM_MAX)) + 
  guides(colour=FALSE) +
  theme_minimal()
plt4

grid.arrange( plt1, plt2)
grid.arrange( plt3, plt4)

