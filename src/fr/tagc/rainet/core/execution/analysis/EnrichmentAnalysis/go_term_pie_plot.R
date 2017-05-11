
## 11-May-2016 Diogo Ribeiro
## Script to a pie plot with GO terms distribution on enriched complexes

library(data.table)
library(plyr)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)

inputFile = "/home/diogo/Documents/RAINET_data/macromolecular_complex_datasets/Wan2015/proportion_wan_go_term_annotations_biological_process_1.txt"

data <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")


data$PERC = round(data$PROPORTION * 100.0,1)

# cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#CC6666", "#9999CC", "#66CC99", "black")
cbPalette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#ffed6f","#ccebc5","#bc80bd")
# cbPalette <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
)


plt1 <- ggplot( data, aes(x = "", y = PROPORTION, fill = factor(GO_NAME)) )  +
  geom_bar(stat = "identity", width = 1, color = "black", size = 0.7, alpha = 0.9) +
  coord_polar(theta="y") +
  # geom_text(aes(y = PERC, label = factor(PERC)), position=position_stack()) +
  scale_fill_manual(name = "GO term", values=cbPalette) +
  blank_theme +
  theme(axis.text.x=element_blank(), legend.position="top", legend.text=element_text(size=13)) +
  guides(fill=guide_legend(nrow=4,byrow=TRUE))
plt1

