# ########################################################################
# This scripts launch the Sweave report that produces statistics
# ########################################################################

library(knitr)

# Get the arguments from the launch command line
args <- commandArgs(TRUE)

# Test if we have enough arguments
if( length(args) != 2){
  stop("Bad argument number")
}

working_dir = args[1]
annotation_input_file = args[2]

knit2pdf('GTEx_statistics.Rnw')
