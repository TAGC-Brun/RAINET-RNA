# ########################################################################
# This scripts launch the Sweave report that produces statistics
# ########################################################################

library(knitr)

# Get the arguments from the launch command line
args <- commandArgs(TRUE)

# Test if we have enough arguments
if( length(args) != 7){
  stop("Rscript: Bad argument number")
}

working_dir = args[1]
output_folder = args[2]
annotation_input_file = args[3]
expression_input_file = args[4]
tx_expression_input_folder = args[5]
tx_expression_avg_file = args[6]
tx_expression_avg_no_outliers_file = args[7]

knit2pdf('GTEx_statistics.Rnw')

file.copy('GTEx_statistics.pdf', output_folder)

