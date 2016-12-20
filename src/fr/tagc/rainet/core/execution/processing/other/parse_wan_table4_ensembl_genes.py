
#20-Dec-2016
# Script to specifically parse Wan et al 2015 http://www.nature.com/nature/journal/v525/n7569/full/nature14877.html supplementary table 4,
# formatting the table for insertion into RAINET DB.
# Using mapping between Ensembl Gene IDs (unique identifiers) instead of Protein names.

workingFolder = "/home/diogo/Documents/RAINET_data/macromolecular_complex_datasets/Wan2015/"

inputFile = workingFolder + "/supp_table4.txt"
conversionFile = workingFolder + "/uniprot-ensembl_gene_id_mapping.tsv"

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
 
with open(inputFile) as inFile:
    header = inFile.readline()
    setComplexes = set()
    for line in inFile:
        line = line.strip()
         
        spl = line.split("\t")
         
        complexID = int( spl[0])
        numberItems = int( spl[1])
        symbols = spl[2] # ensembl gene ids
         
        symbolSpl = symbols.split(";")
         
        assert len(symbolSpl) == numberItems
 
        setComplexes.add( complexID)
         
        for symbol in symbolSpl:
            if symbol in conversionDict:
                found.add( symbol)
                uniprotACs = conversionDict[ symbol]
                for uniprotAC in uniprotACs:
                    outputFile.write( "%s\t%s\n" % (complexID, uniprotAC))
             
    for complexID in sorted( setComplexes):
        outputFileList.write( "%s\n" % (complexID))
 
print "Found mapping for %s IDs" % len( found)

print "FINISHED!"
