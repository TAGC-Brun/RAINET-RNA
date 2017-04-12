
library(data.table)
library(plyr)
require(ggplot2)
require(grid)
require(gridExtra)
library(RColorBrewer)
#source("/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/analysis/RBPDomain/Rscripts/r_functions.R")

# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/Cutoff50/HavugimanaR10000Expr1.0/enrichment_results_filtered.tsv"
# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/Cutoff50/HavugimanaR10000Expr1.0/minProt5/enrichment_results_filtered.tsv"
inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/Cutoff50/HavugimanaR10000Expr1.0/minProt5_topEnrich5/enrichment_results_filtered.tsv"

# inputFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/enrichmentAnalysisStrategy/real/lncRNAs/Cutoff50/NetworkModuleR10000Expr1.0/enrichment_results_filtered.tsv"

data <- fread(inputFile, stringsAsFactors = FALSE, header = TRUE, sep="\t")

###############################
# Data processing
###############################

# calculate maximum number_observed_interactions for each complex
interm = aggregate(data$number_observed_interactions, by = list(data$annotID), max)
colnames(interm) = c("annotID","max_interacting_proteins")
# merge column into original dataset
interm2 <- merge( data, interm,by=c("annotID"))
# add column with number of enrichments
interm2.1 = transform(data, number_enrichments = ave(seq(nrow(data)), annotID, FUN=length))
interm2.2 = subset(interm2.1, select=c("annotID","number_enrichments"))
interm2.3 = unique(interm2.2)
interm2.4 <- merge( interm2, interm2.3,by=c("annotID"))
interm2.5 = subset(interm2.4, select=c("annotID","number_enrichments"))
interm2.6 = unique(interm2.5)
# add number of enrichments (lncRNAs) with max value of number_observed_interactions
interm3 = count(interm2.4, c("annotID","number_observed_interactions","max_interacting_proteins"))
interm3$rnas_with_max = interm3$freq
interm3$freq = NULL
# dataframe with one line per complex
interm4 = interm3[interm3$number_observed_interactions == interm3$max_interacting_proteins, ] 
# select wanted columns from dataframe
interm5 <- merge( data, interm4 ,by=c("annotID"))
interm5 = subset(interm5, select=c("annotID","number_possible_interactions","max_interacting_proteins","rnas_with_max"))
# keep only unique rows
interm6 = unique(interm5)
# merge again number_enrichments column
interm7 = merge( interm6, interm2.6, by=c("annotID"), all.x = TRUE)
# calculate percentage of enrichments where interacting proteins matches the max possible interactions
interm7$perc_max_enrich = round( interm7$rnas_with_max * 100.0 / interm7$number_enrichments, 2)

results = interm7

###############################
# Plotting
###############################

###############################
# Scatter plot of number of maximum protein interactions (in all enrichments of complex) per complex size, with linear regression
###############################
fit = lm(max_interacting_proteins ~ number_possible_interactions, data = results)

### ATTENTION TO XLIM AND YLIM
# Note that position jitter puts dots out of place (can be above limit number of proteins - ab line)
plt1 <- ggplot( results, aes(x = number_possible_interactions, y = max_interacting_proteins) )  +
  scale_color_gradientn(colours = colorRampPalette(c("lightblue", "darkblue", "red"), bias = 1)(5)) +
  scale_size_continuous(range = c(1.5, 7)) + 
  geom_point( aes(color = perc_max_enrich, size = number_enrichments), position="jitter") +  #log10(rnas_with_max))
  geom_text( aes(label = annotID)) +
  geom_abline(intercept = 0, slope = 1) + 
  geom_abline(intercept = signif(fit$coef[[1]],3 ), slope = signif(fit$coef[[2]], 3), color = "brown") + 
  xlim( c(0,50)) +
  ylim( c(0,30)) +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 3),
                     "Intercept =",signif(fit$coef[[1]],3 ),
                     " Slope =",signif(fit$coef[[2]], 3),
                     " P =",signif(summary(fit)$coef[2,4], 3))) +
  theme_minimal()
plt1


# To see if Network modules 'harboring' complexes have similar of different
# results$ratio = results$max_interacting_proteins / results$number_possible_interactions
# 
# results$inside = results$annotID %in% c(4,8,14,33,49,52,56,69,74,87,107,116,178,224)
# 
# ggplot( results, aes(x = ratio, fill = inside)) +
#   geom_histogram( binwidth = 0.01)
# 
# ggplot( results[results$inside == "TRUE"], aes(x = ratio, fill = inside)) +
#   geom_histogram( binwidth = 0.01)
# ggplot( results[results$inside == "FALSE"], aes(x = ratio, fill = inside)) +
#   geom_histogram( binwidth = 0.01)
# 
# ks.test( results[results$inside == "TRUE"]$ratio, results[results$inside == "FALSE"]$ratio)


# ## Scatter plot of number of observed protein interactions (in all enrichments of complex) per complex size, with linear regression
# fit = lm(number_observed_interactions ~ number_possible_interactions, data = data)
# nrow(data)
# 
# plt2 <- ggplot( data, aes(x = number_possible_interactions, y = number_observed_interactions) )  +
#   geom_point( size = 1.5) + 
#   geom_smooth(method=lm, formula = y ~ x, se=FALSE) + 
#   geom_abline(intercept = 0, slope = 1) + 
#   xlim( c(0,50)) +
#   labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 3),
#                      "Intercept =",signif(fit$coef[[1]],3 ),
#                      " Slope =",signif(fit$coef[[2]], 3),
#                      " P =",signif(summary(fit)$coef[2,4], 3)))
# plt2
# 


# ### wanted dataframe: AnnotID\tAnnot size\tEnrichments
# 
# # add column with number of enrichments
# interm = transform(data, number_enrichments = ave(seq(nrow(data)), annotID, FUN=length))
# # select wanted columns from dataframe
# interm2 = subset(interm, select=c("annotID","number_possible_interactions","number_enrichments"))
# # keep only unique rows, since each complex will only have a value for number enrichments and possible interactions
# results = unique(interm2)
# 
# nrow(results)
# 
# ### Plotting
# 
# ## Scatter plot with linear regression
# fit = lm(number_enrichments ~ number_possible_interactions, data = results)
# 
# plt1 <- ggplot( results, aes(x = number_possible_interactions, y = number_enrichments) )  +
#   geom_point( size = 1.5) + 
#   geom_smooth(method=lm, formula = y ~ x, se=FALSE) + 
#  # xlim( c(20,120)) +
#   labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 3),
#                      "Intercept =",signif(fit$coef[[1]],3 ),
#                      " Slope =",signif(fit$coef[[2]], 3),
#                      " P =",signif(summary(fit)$coef[2,4], 3)))
# plt1
# 
# 
# ## Histogram of number of enrichments per number of interacting proteins
# plt2 <- ggplot( data, aes(x = number_observed_interactions) )  +
#   geom_histogram( binwidth = 1) + 
#   xlim(c(2,25))
# plt2
