
# 2017 Diogo Ribeiro
# Script to plot/table results from fisher exact test and its odds ratio on groups of lncRNAs

library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)

##############################
# Basic one-way all vs all table
##############################

# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/outFile.tsv"
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/experimental_interactions/outFile.tsv"
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/rbp_enrichments/outFile.tsv"
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/rbp_enrichments/outFile_complex_dataset.tsv"
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/outFile_simple.tsv"
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/outFile_complex_dataset.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/structure_comparison/outFile.tsv"
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/structure_comparison/outFile_gencodebasic_background.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/minProt_topEnrich/outFile_minProt_topEnrich.tsv"
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/RBP_analysis/GroupOddsRatio/cutoff50_background_rbp/rbp.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

# get all different categories  
categories = unique(dataset$ExternalList)

# create a table for each category
for (i in categories){
  grid.newpage()
  grid.table( dataset[dataset$ExternalList == i])
}

# Whole summary
grid.newpage()
grid.table( dataset)


##############################
# Two-sided colored table (as in Mukherjee2016)
##############################
### change format of dataset: one line per external dataset and 'yes' or 'no', value is odds ratio

# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/outFile_simple.tsv_two_sided.tsv"
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/rbp_enrichments/outFile.tsv_two_sided.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/minProt_topEnrich/outFile_minProt_topEnrich.tsv_two_sided.tsv"
dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

### filter dataset to have only scaffolding candidates
# Proteome-wide
# dataset = dataset[dataset$TranscriptGroup == "2-scaffolding_candidates"]
# RBP-only
# dataset = dataset[dataset$TranscriptGroup == "2-RBP-enriched_lncRNAs"]
# MinProt topEnrich parameter
dataset = dataset[dataset$TranscriptGroup == "minProt5_topEnrich5"]


# # excluding Mukherjee2016 because it has infinite odds ratio
# dataset = dataset[dataset$ExternalList != "Mukherjee2016"]
# excluding Necsulea2014 because there is not significance
#dataset = dataset[dataset$ExternalList != "Necsulea2014"]


## make Inf odds ratio appear as NA, so that we can put a good enrichment color on it
dataset[dataset$OddsRatio == "Inf"]$OddsRatio = "NA"

plt1 = ggplot( dataset, aes(x = ExternalList, y = InGroup)) +
  geom_tile( aes( fill = OddsRatio), colour = "black", size = 1) + 
  scale_fill_continuous( low = "white", high = "#de2d26", name = "Fisher's Exact Test \nOdds ratio", na.value = "#de2d26") +  
  xlab("Orthogonal lncRNA gene dataset") +
  ylab("Scaffolding candidate") +
  geom_text( label = dataset$Overlap, size = 8) +
  theme_minimal() + 
  theme(text = element_text(size=20)) +
  theme(axis.title.x=element_text(vjust=-0.6))
plt1

# print as 20.91 x 4.50 inches

##############################
# Colored table for complex datasets
##############################

inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/LncRNAGroupAnalysis/LncRNAGroupOddsRatio/real/outFile_complex_dataset.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

plt2 = ggplot( dataset, aes(x = TranscriptGroup, y = ExternalList)) +
  geom_tile( aes( fill = OddsRatio), colour = "black", size = 1) +
  scale_fill_continuous( low = "white", high = "#de2d26", name = "Fisher's Exact Test \nOdds ratio", na.value = "#de2d26", limits = c(1.,2.)) +
  scale_x_discrete( labels = c("Corum","Wan 2015","BioPlex","Network modules","All candidates")) +
  xlab("Protein dataset") +
  ylab("Functional/Conserved lncRNAs") +
  geom_text( label = dataset$Overlap, size = 8) +
  theme_minimal() + 
  theme(text = element_text(size=20)) +
  theme(axis.title.x=element_text(vjust=-0.6), axis.text.y=element_blank())
plt2

# print as 20.91 x 4.50 inches

