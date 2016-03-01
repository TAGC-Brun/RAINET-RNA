# ########################################################################
# This scripts launch the Sweave report that produces statistics
# for RAINET Analysis Strategy
# ########################################################################

library(knitr)

# Get the arguments from the launch command line
args <- commandArgs(TRUE)

# Test if we have enough arguments
if( length(args) != 5){
  stop("Rscript: Bad argument number")
}

working_dir = args[1]
sweave_file = args[2]
output_folder = args[3]
parameters_log = args[4]
report_rna_numbers = args[5]

knit2pdf(sweave_file)
