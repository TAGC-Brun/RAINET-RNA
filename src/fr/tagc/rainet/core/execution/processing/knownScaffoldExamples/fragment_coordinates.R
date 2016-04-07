
library(ggplot2)
library(gplots)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)
library(data.table)

inputFile = "/home/diogo/Documents/RAINET_data/catRAPID/webserver_results/XIST_MOUSE/Q4VBD9_and_Q99JB6_fragments.txt"

dataset = read.table( inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

CATRAPID_ZS_MEAN = 23.25
CATRAPID_ZS_STD = 37.90

plt1 <- ggplot(dataset) +
  labs(x="XIST transcript coordinate",y="score of fragment") +
  theme_minimal()
#for each fragment
for (i in seq(1:length(dataset$ini))){
  ini = dataset$ini[i]
  end = dataset$end[i] 
  score = dataset$score[i] * CATRAPID_ZS_STD - CATRAPID_ZS_MEAN
  prot = dataset$prot[i]
  if (prot == "Q4VBD9"){
    cor = "'red'"
  }else{
    cor = "'blue'"
  }
  loop_input = paste("annotate('segment', x = ",ini,", xend = ",end,", y = ",score,", yend =",score,", colour =",cor,", size = 2)" )
  plt1 <- plt1 + eval(parse(text=loop_input))  
  }

# # RepA
# plt1 <- plt1 + annotate('segment', x = 227, xend = 760, y = max( dataset$score) + 0.5, yend = max( dataset$score) + 0.5, colour = "green", size = 6)
# # 4R
# plt1 <- plt1 + annotate('segment', x = 3098, xend = 4713, y = max( dataset$score) + 0.5, yend = max( dataset$score) + 0.5, colour = "green", size = 6)
# 
plt1
