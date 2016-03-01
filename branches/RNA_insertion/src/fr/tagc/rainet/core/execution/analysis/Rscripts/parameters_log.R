
#inputFile = "/home/diogo/workspace/tagc-rainet-RNA/test/fr/tagc/rainet/core/test_results/parameters.log"
inputFile = paste( output_folder, parameters_log, sep = "/")

inputFile

inputData <- read.table(inputFile, header = TRUE, sep = "\t")

print( inputData)
