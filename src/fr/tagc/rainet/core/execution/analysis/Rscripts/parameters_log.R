
require(grid)
require(gridExtra)

#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/parameters.log"
inputFile = paste( output_folder, parameters_log, sep = "/")

inputData <- read.table(inputFile, header = TRUE, sep = "\t")

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.7)),
  colhead = list(fg_params=list(cex = 0.8)),
  rowhead = list(fg_params=list(cex = 0.8)))

myt <- gridExtra::tableGrob( inputData, theme = mytheme)

grid.draw(myt)
