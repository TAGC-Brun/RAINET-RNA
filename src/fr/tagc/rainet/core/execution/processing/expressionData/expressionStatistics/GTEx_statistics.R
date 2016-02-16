# ########################################################################
# This scripts launch the Sweave report that produces statistics
# ########################################################################

library(knitr)

# Get the arguments from the launch command line
args <- commandArgs(TRUE)

# Test if we have enough arguments
if( length(args) != 5){
  stop("Bad argument number")
}

working_dir = args[1]
annotation_input_file = args[2]
expression_input_file = args[3]
expression_sample_input_file = args[4]
tx_expression_input_folder = args[5]

knit2pdf('GTEx_statistics.Rnw')
