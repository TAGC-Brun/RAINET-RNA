
# 22-Mar-2017 Diogo Ribeiro
# Script to plot a venn diagram comparing lncRNA scaffold candidates between proteome-wide and RBP-only enrichment analysis approaches

library(data.table)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
library(VennDiagram)
#source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

##########################################
# Parameters
##########################################
proteomeApproachFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/Cutoff50/merge/list_proteome_enriched_rnas_havugimana2012_bioplex_wan2015_networkModules.txt"
rbpApproachFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/Cutoff50_rbps/merge/list_rbp_enriched_rnas_havugimana2012_bioplex_wan2015_networkModules.txt"
background = 15230 # Background is for example the 15230 tested transcripts

proteomeApproachList = fread(proteomeApproachFile, stringsAsFactors = FALSE, header = FALSE)
rbpApproachList = fread(rbpApproachFile, stringsAsFactors = FALSE, header = FALSE)

# Convert lists to set
proteomeApproachSet = c(proteomeApproachList$V1)
rbpApproachSet = c(rbpApproachList$V1)

##########################################
# basic set operations
##########################################
setUnion = union( proteomeApproachSet, rbpApproachSet)
setIntersect = intersect( proteomeApproachSet, rbpApproachSet)
setProteomeMinusRBP = setdiff( proteomeApproachSet, rbpApproachSet)
setRBPMinusProteome = setdiff( rbpApproachSet, proteomeApproachSet)

paste( "Proteome dataset length:", length(proteomeApproachSet))
paste( "RBP dataset length:", length(rbpApproachSet))
paste( "Set union length:", length(setUnion))
paste( "Set intersection length:", length(setIntersect))
paste( "Set Proteome minus RBP length:", length(setProteomeMinusRBP))
paste( "Set RBP minus Proteome length:", length(setRBPMinusProteome))


##########################################
# Draw venn diagram between two sets
##########################################
grid.newpage()
# Arguments: area1, area2, cross.area. The individual areas will be substrated by the cross area
draw.pairwise.venn( length(proteomeApproachSet), length(rbpApproachSet), length(setIntersect),
     category = c("Proteome candidates", "RBP candidates"), 
     # lty = rep("blank", 2), 
     fill = c("#d8b365", "#5ab4ac"),
     alpha = rep(0.7, 2), 
     cat.pos = c(-15, 15), 
     cex = rep(2,3),
     cat.cex = rep(1.5, 2),
     # cat.dist = rep(0.025, 2)
  )


##########################################
# Statistical significance of intersection.
##########################################

## Example contigency table
#            | In proteome | not in proteome 
#  in RBP    |     572     |      216
# not in RBP |     281     |   15230-1069

ft = fisher.test( matrix(c( length(setIntersect), length(setRBPMinusProteome), length(setProteomeMinusRBP), background - length( union)), nrow = 2))

grid.text(paste("Intersection\nFisher's exact test:\np-val: ", ft$p.value, "\nOdds ratio: ", round(ft$estimate,2)), x = unit(0.5, "npc"), y = unit(0.7, "npc"))

