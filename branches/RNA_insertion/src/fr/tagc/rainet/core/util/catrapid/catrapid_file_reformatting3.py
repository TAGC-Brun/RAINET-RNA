
# 2-May-2017 Diogo Ribeiro

# Simple routine to reformat from:
# sp|Q53F19|NCBP3_HUMAN ENST00000423712   -58.53  0.00    0.00
# To:
# ENST00000423712    Q53F19    -58.53

import sys

inputFile = sys.argv[1]
outputFile = open( sys.argv[2], "w")

with open( inputFile, "r") as inFile:
    for line in inFile:
        spl = line.strip().split( "\t")
        
        ids = spl[0]
        prot, txID = ids.split(" ")
        score = spl[1]
        otherColums = spl[1:]

        protID = prot.split("|")[1]
        
        outputFile.write("%s\t%s\t%s\n" % (txID, protID, "\t".join(otherColums)) )
        
outputFile.close()
print "FINISHED!"
