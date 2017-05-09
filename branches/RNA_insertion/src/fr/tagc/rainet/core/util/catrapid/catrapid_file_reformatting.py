
# 13-Feb-2017 Diogo Ribeiro

# Simple routine to reformat from:
# ENST00000448096|P52298  90.87   1
# To:
# sp|O00425|sp ENST00000470491    1

import sys

inputFile = sys.argv[1]
outputFile = open( sys.argv[2], "w")

with open( inputFile, "r") as inFile:
    for line in inFile:
        spl = line.strip().split( "\t")
        
        txID,protID = spl[0].split("|")
        score = spl[1]
        otherColums = spl[1:]
        
        outputFile.write("sp|%s|sp %s\t%s\n" % (protID, txID, "\t".join(otherColums) ) )
        
outputFile.close()
print "FINISHED!"