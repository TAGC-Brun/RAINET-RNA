# ########################################################################
# This scripts launch the Sweave report that produces statistics
# for RAINET Analysis Strategy
# ########################################################################

library(knitr)

# Get the arguments from the launch command line
args <- commandArgs(TRUE)

# Test if we have enough arguments
if( length(args) != 8){
  stop("Rscript: Bad argument number")
}

working_dir = args[1]
sweave_file = args[2]
output_folder = args[3]
parameters_log = args[4]
report_rna_numbers = args[5]
report_interaction_numbers = args[6]
report_interaction_scores = args[7]
report_interaction_partners = args[8]

knit2pdf(sweave_file)
