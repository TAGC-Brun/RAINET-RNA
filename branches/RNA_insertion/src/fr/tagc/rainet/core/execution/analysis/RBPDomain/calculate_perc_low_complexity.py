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

outFile.write( "uniprotac\ttotal_length\tlowercase\tuppercase\ttag\tperc_low\n")

# for each entry, calculate percentage underscore, write to file
for ID in fastaData:
    seq = fastaData[ ID]
    lowerN = n_lower_chars( seq)
    upperN = n_upper_chars( seq)
    perc = lowerN * 100.0 / (lowerN + upperN) 
    
    # measure presence of low complexity in the terminus region (first 25 and last 25 aa)

    # From Coletta et al. 2010
    # t-LCRs were defined as regions starting or ending at no more than 25 amino acids from either sequence
    # extremity, c-LCRs as regions starting or ending at least 50 amino acids from either sequence extremity.
    
    terminus = seq[0:25] + seq[-25:]
    central = seq[50:-50]

    assert len(terminus) == 50
    assert len(central) <= max(len(seq) - 100, 0)

    # location tags:
    # 0 - no LC
    # 1 - terminus LC
    # 2 - central LC
    # 3 - other LC

    terminusBool = n_lower_chars( terminus) > 0
    centralBool = n_lower_chars( central) > 0
    otherBool = n_lower_chars( seq) > 0

    tag = -1
    if terminusBool and centralBool:
        tag = "other LC"
    elif terminusBool == 0 and centralBool == 0 and otherBool:
        tag = "other LC"
    elif terminusBool and centralBool == 0:
        tag = "terminus LC"
    elif terminusBool == 0 and centralBool:
        tag = "central LC"
    else:
        tag = "no LC"

    assert lowerN + upperN == len(seq), "number of lower and upper characters should sum up to protein length, if not there may be wrong characters."
    assert perc >= 0 and perc <= 100
    
    outFile.write ("%s\t%i\t%i\t%i\t%s\t%.2f\n" % (ID, len(seq), lowerN, upperN, tag, perc) )

outFile.close()

print "FINISHED!"



