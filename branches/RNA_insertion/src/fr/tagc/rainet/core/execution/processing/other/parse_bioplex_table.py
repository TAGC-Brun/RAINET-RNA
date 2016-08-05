
#05-Aug-2016
# Script to specifically parse Bioplex complexes dataset
# formatting the table for insertion into RAINET DB.

workingFolder = "/home/diogo/Documents/RAINET_data/macromolecular_complex_datasets/Bioplex"

inputFile = workingFolder + "/bioplex_complexes.txt"
conversionFile = workingFolder + "/uniprot-yourlist%3AM2016080583C3DD8CE55183C76102DC5D3A26728BA8BDD02.tab"

outputFile = open( workingFolder + "/bioplex_complexes_annotation.txt", "w" )
outputFileList = open( workingFolder + "/bioplex_complexes_list.txt", "w" )

### read conversion file

conversionDict = {} # key -> ID from Wan, val -> list of uniprotAC

with open( conversionFile) as inFile:
    header = inFile.readline()
    for line in inFile:
        line = line.strip()
        
        spl = line.split("\t")
        
        uniprotAC = spl[0]
        
        otherIDs = spl[1].split(" ")
                
        for ID in otherIDs:
            if ID not in conversionDict:
                conversionDict[ ID] = []
            conversionDict[ ID].append( uniprotAC)

### read annotation file
 
found = set()
 
with open(inputFile) as inFile:
    header = inFile.readline()
    setComplexes = set()
    for line in inFile:
        line = line.strip()
         
        spl = line.split("\t")
         
        complexID = spl[0]
        symbol = spl[2]
          
        setComplexes.add( complexID)
         
        if symbol in conversionDict:
            found.add( symbol)
            uniprotACs = conversionDict[ symbol]
            for uniprotAC in uniprotACs:
                outputFile.write( "%s\t%s\n" % (complexID, uniprotAC))
             
    for complexID in sorted( setComplexes):
        outputFileList.write( "%s\n" % (complexID))
 
print "Found mapping for %s IDs" % len( found)

print "FINISHED!"
