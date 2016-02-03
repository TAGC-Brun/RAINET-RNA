# ########################################################################################################################
# This script aims to load a PPI network file (two column file, one protein per column) and remove from it the self loop
# i.e. the interaction that goes from a protein to itself.
# The script takes as first argument the path to the input PPI network file and as second argument the path to the
# desired output file (same format as input file)
# ########################################################################################################################

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
network_input_path=as.character( args[1])
network_output_path=as.character( args[2])

# Read the network (two columns file)
network_df = read.table( file = network_input_path, stringsAsFactors = FALSE, sep="\t", header = FALSE)

# Locate the lines with self loop
self_loop_indexes = which( network_df$V1 == network_df$V2)

# Build the network df with no loop
network_noself_df = network_df[ -self_loop_indexes,]

# Export result to fil
write.table( network_noself_df, file = network_output_path, sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
