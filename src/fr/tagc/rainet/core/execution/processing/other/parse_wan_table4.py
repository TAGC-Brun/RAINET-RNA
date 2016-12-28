
#05-Aug-2016
# Script to specifically parse Wan et al 2015 http://www.nature.com/nature/journal/v525/n7569/full/nature14877.html supplementary table 4,
# formatting the table for insertion into RAINET DB.

workingFolder = "/home/diogo/Documents/RAINET_data/macromolecular_complex_datasets/Wan2015/"

inputFile = workingFolder + "/supp_table4.txt"
conversionFile = workingFolder + "/uniprot-yourlist%3AM20160805325D09DDFD8B5D0CFB8A8926E064CD7EB15F2BQ.tab"

outputFile = open( workingFolder + "/Wan2015_s4_complexes_annotation.txt", "w" )
outputFileList = open( workingFolder + "/Wan2015_s4_complexes_list.txt", "w" )

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

duplicates = set()
 
with open(inputFile) as inFile:
    header = inFile.readline()
    setComplexes = set()
    for line in inFile:
        line = line.strip()
         
        spl = line.split("\t")
         
        complexID = int( spl[0])
        numberItems = int( spl[1])
        ensgs = spl[2]
        symbols = spl[3]
         
        assert len( ensgs.split(";")) == numberItems
         
        symbolSpl = symbols.split(";")
         
        assert len(symbolSpl) == numberItems
 
        setComplexes.add( complexID)
         
        for symbol in symbolSpl:
            if symbol in conversionDict:
                found.add( symbol)
                uniprotACs = conversionDict[ symbol]
                
                # for each corresponding uniprotAC.. note that this can produce duplicate mappings
                for uniprotAC in uniprotACs:
                    outputFile.write( "%s\t%s\n" % (complexID, uniprotAC))

                if len( uniprotACs) > 1:
                    print "DUPLICATE:", complexID, symbol, uniprotACs
                    duplicates.add( symbol)

             
    for complexID in sorted( setComplexes):
        outputFileList.write( "%s\n" % (complexID))
 
print "Found mapping for %s IDs" % len( found)
print "Duplicate ID mapping for %s IDs" % len( duplicates)

print "FINISHED!"
