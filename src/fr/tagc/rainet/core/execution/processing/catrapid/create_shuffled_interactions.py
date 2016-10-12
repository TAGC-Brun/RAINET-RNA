
import sys
import random

#12-Oct-2016
# Simple script to create set of catRAPID interactions, by replacing/shuffling names of proteins.
# The number and values of scores for an RNA will be kept, only their binding partners change.

proteinListFile = sys.argv[1]
catrapidFile = sys.argv[2]

### read list of proteins and swap names of proteins

proteinList = [line.strip() for line in open( proteinListFile,"r").readlines()]

randomizedProteinList = random.sample( proteinList, len( proteinList))

assert( len( set(proteinList)) == len( set(randomizedProteinList)) )

swapper = {}

for i in xrange( 0, len( proteinList)):
    swapper[ proteinList[ i]] = randomizedProteinList[ i]

assert ( len(swapper) == len( proteinList) )

### read catrapid file and modify the protein names
 
outFile = open( catrapidFile + "_shuffled","w")

count = 0

with open( catrapidFile,"r") as inFile:
    for line in inFile:
        
        count+=1
        if count % 10000000 == 0:
            print "Processed %s lines.." % count
        
        spl = line.split(" ")
        protID = spl[0]
        newProtID = swapper[ protID]
        
        newLine = newProtID + " " + spl[1]

        outFile.write( newLine) 
 
outFile.close()

print "FINISHED!"