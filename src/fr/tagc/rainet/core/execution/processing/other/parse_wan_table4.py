
#05-Aug-2016
# Script to specifically parse Wan et al 2015 http://www.nature.com/nature/journal/v525/n7569/full/nature14877.html supplementary table 4,
# formatting the table for insertion into RAINET DB.

workingFolder = "/home/diogo/Documents/RAINET_data/macromolecular_complex_datasets/Wan2015/"

inputFile = workingFolder + "/supp_table4.txt"
conversionFile = workingFolder + "/uniprot-yourlist%3AM20160805325D09DDFD8B5D0CFB8A8926E064CD7EB15F2BQ.tab"
conversionFileGene = workingFolder + "/uniprot-ensembl_gene_id_mapping.tsv"

outputFile = open( workingFolder + "/Wan2015_s4_complexes_annotation.txt", "w" )
outputFileList = open( workingFolder + "/Wan2015_s4_complexes_list.txt", "w" )

###


### read conversion file between ENSG and uniprotAC

conversionDictGene = {} # key -> ID from Wan, val -> list of uniprotAC

with open( conversionFileGene) as inFile:
    header = inFile.readline()
    for line in inFile:
        line = line.strip()
        
        spl = line.split("\t")
        
        uniprotAC = spl[0]
        
        otherIDs = spl[1].split(" ")
                
        for ID in otherIDs:
            if ID not in conversionDictGene:
                conversionDictGene[ ID] = []
            conversionDictGene[ ID].append( uniprotAC)


### read conversion file between Protein name and uniprotAC

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
## ID mapping approach
# Map using protein names, as they map to a large number of uniprotACs.
# When there is same protein name mapping to several uniprotACs, try to use ENSG mapping to make them unique
 
found = set()
complexProteinPair = set() # store complex-protein pairs to write to file, avoiding duplicate entries

with open(inputFile) as inFile:
    header = inFile.readline()
    setComplexes = set()
    for line in inFile:
        line = line.strip()
         
        spl = line.split("\t")
         
        complexID = int( spl[0])
        numberItems = int( spl[1])
        ensgs = spl[2]
        ensgSpl =  ensgs.split(";")
         
        assert len( ensgSpl) == numberItems

        # get list of uniprotACs found with the ENSG
        genePool = set()
        for ensg in ensgSpl:
            if ensg in conversionDictGene:
                uniIDs = conversionDictGene[ ensg]
                for uni in uniIDs:
                    genePool.add( uni)

        symbols = spl[3]         
        symbolSpl = symbols.split(";")
         
        assert len(symbolSpl) == numberItems
 
        setComplexes.add( complexID)

        symbolsPool = set()         
        for symbol in symbolSpl:
            if symbol in conversionDict:
                found.add( symbol)
                uniprotACs = conversionDict[ symbol]
                
                added = 0
                # if there are duplicate mappings, include only the ones also present with the ENSG ID
                if len( uniprotACs) > 1:
                    for uniprotAC in uniprotACs:
                        if uniprotAC in genePool:
                            added +=1
                            tag = str(complexID) + "|" + uniprotAC
                            complexProteinPair.add( tag)

                    if added > 1:
                        print "Duplicate in both protein name and ENSG", complexID, symbol, uniprotACs

                else:
                    # there is only one, loop just because item is a list
                    for uniprotAC in uniprotACs:
                        tag = str(complexID) + "|" + uniprotAC
                        complexProteinPair.add( tag)

    # Write to file
    for complexPair in sorted( complexProteinPair):
        complexPairSpl = complexPair.split("|")
        outputFile.write( "%s\t%s\n" % (complexPairSpl[0], complexPairSpl[1]))   
             
    for complexID in sorted( setComplexes):
        outputFileList.write( "%s\n" % (complexID))
 
print "Found mapping for %s IDs" % len( found)

print "FINISHED!"
