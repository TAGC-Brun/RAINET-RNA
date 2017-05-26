
## 11-May-2016 Diogo Ribeiro
## Script to create a pie plot with GO terms distribution on enriched complexes

library(data.table)
library(plyr)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)

# Wan 2015
# inputFile = "/home/diogo/Documents/RAINET_data/macromolecular_complex_datasets/Wan2015/proportion_wan_go_term_annotations_biological_process_1.txt"

# Network Modules
inputFile = "/home/diogo/Documents/RAINET_data/macromolecular_complex_datasets/NetworkModules/proportion_wan_go_term_annotations_biological_process_1.txt"

data <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")


data$PERC = round(data$PROPORTION * 100.0,1)

# # Colors for Wan 2015
# cbPalette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#ffed6f","#ccebc5","#bc80bd")

# Colors for Network modules (to match Wan 2015 when possible)
cbPalette <- c("#ffffb3","#bebada","#fb8072","#80b1d3", "#4eb3d3","#fdb462","#b3de69","#d9d9d9","#ffed6f","#ccebc5","#bc80bd")


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
)

plt1 <- ggplot( data, aes(x = "", y = PROPORTION, fill = GO_NAME) )  +
  geom_bar( width = 1, color = "black", size = 0.7, alpha = 0.9) + #stat = "identity"
  coord_polar(theta="y") +
  # geom_text(aes(y = PERC, label = factor(PERC)), position=position_stack()) +
  scale_fill_manual(name = "GO term", values=cbPalette) +
  blank_theme +
  theme(axis.text.x=element_blank(), legend.position="top", legend.text=element_text(size=10), legend.key.size = unit(1,"line")) +
  guides(fill=guide_legend(nrow=4,byrow=TRUE))
plt1

#print as 7.00 x 5.91 inches
