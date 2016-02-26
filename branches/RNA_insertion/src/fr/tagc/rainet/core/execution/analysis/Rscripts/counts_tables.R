
library(ggplot2)
library(reshape)
library(RColorBrewer)

inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/rna_numbers.tsv"
# inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/interaction_numbers.tsv"
# inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/rna_expression.tsv"
# inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/Report/rna_expression_data_presence.tsv"

rna_numbers <- read.table(inputFile, header = TRUE, sep = "\t")
rna_numbers <- melt(rna_numbers)

# Re-ordering labels
rna_numbers$Data <- factor(rna_numbers$Data, levels = rev(levels(rna_numbers$Data)))
rna_numbers$variable <- factor(rna_numbers$variable, levels = rev(levels(rna_numbers$variable)))

# Bar plot for numbers of items before and after filter
ggplot( rna_numbers, aes(x = variable, y = value, fill = Data)) + 
  geom_bar( stat = "identity", position=position_dodge()) +
  coord_flip() + 
  xlab( "Class of RNA") +
  ylab( "Frequency")

