
#05-Aug-2016
# Script to specifically parse Wan et al 2015 http://www.nature.com/nature/journal/v525/n7569/full/nature14877.html supplementary table 4,
# formatting the table for insertion into RAINET DB.

workingFolder = "/home/diogo/Documents/RAINET_data/macromolecular_complex_datasets/Wan2015/"

inputFile = workingFolder + "/supp_table4.txt"

outputFile = open( workingFolder + "/Wan2015_s4_complexes_annotation.txt", "w" )
outputFileList = open( workingFolder + "/Wan2015_s4_complexes_list.txt", "w" )

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
            outputFile.write( "%s\t%s\n" % (complexID, symbol))
            
    for complex in sorted( setComplexes):
        outputFileList.write( "%s\n" % (complex))
        
print "FINISHED!"
