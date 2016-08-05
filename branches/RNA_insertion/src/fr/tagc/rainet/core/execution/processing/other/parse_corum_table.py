
#05-Aug-2016
# Script to specifically parse Corum database allComplexesCore.csv file
# formatting the table for insertion into RAINET DB.

workingFolder = "/home/diogo/Documents/RAINET_data/macromolecular_complex_datasets/Corum"

inputFile = workingFolder + "/allComplexesCore.csv"

outputFile = open( workingFolder + "/corum_complexes_annotation.txt", "w" )
outputFileList = open( workingFolder + "/corum_complexes_list.txt", "w" )

WANTED_ORGANISM = "Human"

#### NOTES ####
# from Corum: SwissProt primary accessions of the protein complex subunits. If subunits can not be identified unam biguously, the respective groups of SwissProt entries are presented within brackets. An example is comp$
# I will only retain the ones not in brackets
###############

### read annotation file
 
with open(inputFile) as inFile:
    header = inFile.readline()
    for line in inFile:
        line = line.strip()
         
        spl = line.split(";")
        
        complexID = spl[0]
        complexName = spl[1]
        organism = spl[3]
        uniprotACs = spl[4]
        method = spl[6]

        if organism != WANTED_ORGANISM:
            continue

#        print complexID, complexName, organism, uniprotACs
        
        # write up complex info
        outputFileList.write( "%s\t%s\t%s\n" % (complexID, complexName, method))

        # write up annotation info
        uniprotACsList = uniprotACs.split(",")
        
        boo = 0 # flag for controlling start of parenthesis
        for uniprotAC in uniprotACsList:
            if uniprotAC.startswith("("):
                boo = 1
                continue
            elif boo:
                continue
            else:        
                outputFile.write( "%s\t%s\n" % (complexID, uniprotAC))
             

print "FINISHED!"
