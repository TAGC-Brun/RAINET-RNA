
library(ggplot2)
library(gplots)
library(reshape)
library(RColorBrewer)
require(grid)
require(gridExtra)
library(data.table)

##################################################
### Distribution number proteins per pathway
##################################################

#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/prot_per_annotation.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/prot_per_annotation.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/test/prot_per_annotation.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/CorumR1000Expr1.0/prot_per_annotation.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/non_redundant/HavugimanaR1000Expr1.0/prot_per_annotation.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

paste( "Total number of annotations:" , nrow( dataset))
#paste( "Total number of proteins with annotation:" , sum(dataset$total_prot_in_pathway) )
#paste( "Total number of proteins with interaction data with annotations:" , sum(dataset$total_prot_with_interaction_data) )

#dataset[order (dataset$total_prot_in_pathway, decreasing = TRUE),]

plt0 <- ggplot( dataset) + 
  geom_histogram( binwidth = 1,  aes( x = total_prot_in_pathway), fill = "red", alpha = 0.5, position="identity") + 
  theme_minimal() +
#  xlim( c(0,50)) +
  xlab( "Number of proteins in annotation") +
  ylab( "Annotation frequency")
plt0

plt1 <- ggplot( dataset) + 
  geom_histogram( binwidth = 1,  aes( x = total_prot_with_interaction_data), fill = "blue", alpha = 0.5, position="identity") + 
  theme_minimal() +
#  xlim( c(0,50)) +
  xlab( "Number of interacting proteins in annotation") +
  ylab( "Annotation Frequency")
plt1

# plt1 <- ggplot( dataset) + 
#   geom_density( aes( x = total_prot_with_interaction_data), fill = "blue", alpha = 0.5, position="identity") + 
#   geom_density( aes( x = total_prot_in_pathway), fill = "red", alpha = 0.5, position="identity") + 
# #  scale_fill_manual(values = c("total_prot_with_interaction_data", "total_prot_in_pathway")) + 
#   theme_minimal() +
#   xlab( "Number of interacting proteins in annotation") +
#   ylab( "Annotation Frequency")
# plt1

grid.arrange(plt0, plt1)

##################################################
### Distribution pathways per protein
##################################################

# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/annotation_per_prot.tsv"
#inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/CorumR1000Expr1.0/annotation_per_prot.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/non_redundant/HavugimanaR1000Expr1.0/annotation_per_prot.tsv"

dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

paste( "Total number of proteins with annotation:" , nrow( dataset))
paste( "Total number of proteins with annotation and with interaction:" , sum(dataset$with_interaction_data) )

#dataset[order (dataset$total_annotations, decreasing = TRUE),]
#datasetInter = dataset[ dataset$with_interaction_data == 1]

plt2 <- ggplot( dataset) + 
  geom_histogram( binwidth = 1,  aes( x = total_annotations), fill = "dark green", alpha = 0.5, position="identity") + 
#  geom_histogram(data = datasetInter, binwidth = 1,  aes( x = total_annotations), fill = "blue", alpha = 0.5, position="identity") + 
  theme_minimal() +
  xlab( "Number annotations per protein") +
  ylab( "Protein frequency")
plt2

#dataset[order(dataset$total_annotations)]

# ##################################################
# ### TODO (network modules check)
# ##################################################
# 
# inputFile = "/home/diogo/testing/networkModules/networkModules.tsv"
# 
# dataset <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")
# 
# head( dataset[order( dataset$number_interactions, decreasing = TRUE)], 10)
# 
# head( dataset[order( dataset$number_possible_interactions, decreasing = TRUE)], 10)
