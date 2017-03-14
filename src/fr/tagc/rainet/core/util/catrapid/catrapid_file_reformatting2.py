
# 14-Mars-2017 Diogo Ribeiro

# Simple routine to reformat from:
# Q9Y6A4  ENST00000563770 0.0 (all tab-separated) (this is format of interaction expression data)
# To:
# sp|O00425|sp ENST00000470491    1 (space-separated, then tab-separated)

import sys

inputFile = sys.argv[1]
outputFile = open( sys.argv[2], "w")

count = 0

with open( inputFile, "r") as inFile:
    for line in inFile:
        protID,txID,score = line.strip().split( "\t")

        count +=1
        if count % 10000000 == 0:
            print "Processed %s lines.." % count
        
        outputFile.write("sp|%s|sp %s\t%s\n" % (protID, txID, score) )
        
outputFile.close()
print "FINISHED!"