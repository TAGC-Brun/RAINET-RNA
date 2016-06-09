#!/usr/bin/env python2

# 08-June-2016
# Simple script to calculate percentage of low complexity region in sequences already masked by segmasker

import sys

if len(sys.argv) != 3:
        print ("Usage: python2 %s masked_fasta_file output_file" % sys.argv[0])
        sys.exit(1)

inFile = open( sys.argv[1], "r")
outFile = open( sys.argv[2], "w")

def n_lower_chars(string):
    return sum(1 for c in string if c.islower())

def n_upper_chars(string):
    return sum(1 for c in string if c.isupper())


#read fasta file
fastaData = {}
text = ""
boo = 0
ID = ""
for line in inFile:
    line = line.strip()
    if line.startswith(">"):
        if boo:
            fastaData[ID] = text
        text = ""
        ID = line
        # process ID to have uniprot ac
        # i.e.: >sp|O60543|CIDEA_HUMAN to O60543
        ID = ID.split( "|")[1]
    else:
        boo = 1
        text += line

fastaData[ID] = text #for the last iteration

outFile.write( "uniprotac\ttotal_length\tlowercase\tuppercase\tperc_lower\n")

# for each entry, calculate percentage underscore, write to file
for ID in fastaData:
    seq = fastaData[ ID]
    lowerN = n_lower_chars( seq)
    upperN = n_upper_chars( seq)
    perc = lowerN * 100.0 / (lowerN + upperN) 

    assert lowerN + upperN == len(seq), "number of lower and upper characters should sum up to protein length, if not there may be wrong characters."
    assert perc >= 0 and perc <= 100
    
    outFile.write ("%s\t%i\t%i\t%i\t%.2f\n" % (ID, len(seq), lowerN, upperN, perc) )

outFile.close()

print "FINISHED!"



