# ########################################################################
# This scripts launch the Sweave report that produces statistics
# ########################################################################

library(knitr)

# Get the arguments from the launch command line
args <- commandArgs(TRUE)

# Test if we have enough arguments
if( length(args) != 6){
  stop("Rscript: Bad argument number")
}

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

working_dir = args[1]
annotation_input_file = args[2]
expression_input_file = args[3]
expression_sample_input_file = args[4]
tx_expression_input_folder = args[5]
tx_expression_avg_file = args[6]

knit2pdf('GTEx_statistics.Rnw')
