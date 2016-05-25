
import sys

#25-May-2016
#Simple script to create subset of a catRAPID library, keeping only the entries matching the IDs provided

listFile = sys.argv[1]
libraryFile = sys.argv[2]
proteinBool = int( sys.argv[3]) # 1 if library os protein, 0 if RNA

itemList = {line.strip() for line in open( listFile,"r").readlines()}

outFile = open( libraryFile + "_pruned","w")

#store proteins that are found
found = set()

with open( libraryFile,"r") as inFile:
    for line in inFile:
        if proteinBool:
            ID = line.split("|")[1]
            if ID in itemList:
                outFile.write(line)
                found.add( ID)
        else:
            ID = line.split("\t")[0]
            if ID in itemList:
                outFile.write(line)
                found.add( ID)           

print "Items found", len(found)
notFound = itemList - found
print "Items not found" , len(notFound)

# reasons for not found, can be because they are not reviewd proteins, or > 750 aa etc.
print notFound

outFile.close()