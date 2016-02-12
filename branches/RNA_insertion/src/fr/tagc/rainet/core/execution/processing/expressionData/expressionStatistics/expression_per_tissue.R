
library(ggplot2)

expression_input_file = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/dataFiles/expression.csv"

# Load the data
expression_df = read.table( expression_input_file, stringsAsFactors = FALSE, header = FALSE, sep=",")

