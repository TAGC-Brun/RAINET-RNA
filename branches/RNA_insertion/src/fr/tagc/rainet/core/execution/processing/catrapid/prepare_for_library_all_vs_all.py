import sys
import os
import shutil
import argparse
import glob
import numpy as np

#===============================================================================
# Started 13-May-2016 
# Diogo Ribeiro
# Script to prepare files for running catRAPID all vs all library
DESC_COMMENT = "Script to prepare files for running catRAPID all vs all library\n \
                This will split both RNA and Protein library files and run all combinations between them"
#===============================================================================

#===============================================================================
# General plan:
# 1) Input is RNA and Protein libraries
# 2) Breakdown subsets of RNAs and Proteins into numbered files
# 3) Confirm if every combination is accounted for
#===============================================================================

# constants
SPLIT_FILES_FOLDER = "/split_files"
PREPARE_FOLDER = "/prepared_folders"


# # 
# Split library files into many files based on batch size
def splitFile( library_file, batch_size, output_folder):

    outPath = output_folder + SPLIT_FILES_FOLDER
    
    if not os.path.exists( outPath):
        os.mkdir( outPath)
        
    # remember -> 10 lines for each RNA/protein
    numMolecules = batch_size * 10
    
    inFileName = os.path.basename( library_file)
    
    # read input file
    with open( library_file, "r") as inFile:
        text = ""
        totalLines = 0 # for statistical purposes
        lineCount = 0
        fileCount = 0
        for line in inFile:
            totalLines +=1
            lineCount += 1
            
            text += line
            
            if lineCount == numMolecules:
                with open( outPath + "/" + inFileName + "_" + str( fileCount) + ".lib", "w") as outFile:
                    outFile.write( text)
                fileCount+=1
                lineCount = 0
                text = ""
                

        # write the last batch
        with open( outPath + "/" + inFileName + "_" + str( fileCount) + ".lib", "w") as outFile:
            outFile.write( text)
    
        
        # check if written correct number of files
        expectedNumberFiles = int( totalLines / float( numMolecules)) + 1        
        writtenNumberFiles = len( glob.glob( outPath + "/" + inFileName + "*"))
        assert expectedNumberFiles == writtenNumberFiles

        print "TOTAL LINES:", library_file, totalLines

# #
# Create a folder ready to be run for each combination of protein and RNA split files
def prepareFolders( rna_library, prot_library, output_folder, catrapid_template):
    
    outPath = output_folder + PREPARE_FOLDER
    
    if not os.path.exists( outPath):
        os.mkdir( outPath)
    
    protFileName = os.path.basename( prot_library)
    rnaFileName = os.path.basename( rna_library)
  
    rnaFiles = glob.glob( output_folder + SPLIT_FILES_FOLDER + "/" + rnaFileName + "*")
    protFiles = glob.glob( output_folder + SPLIT_FILES_FOLDER + "/" + protFileName + "*")

    folderCount = 0
    for protFile in protFiles:
        for rnaFile in rnaFiles:
            folderCount += 1
            newFolder = outPath + "/run" + str( folderCount)
            if not os.path.exists( newFolder):
                os.mkdir( newFolder)
            
            # copy template to newFolder
            os.system( "cp -r %s/* %s" % ( catrapid_template, newFolder) )
            
            # copy protein file and rna file to that folder
            os.system( "cp -r %s %s/prot/prot.lib" % ( protFile, newFolder) )
            os.system( "cp -r %s %s/rna/rna.lib" % ( rnaFile, newFolder) )

            if folderCount % 100 == 0:
                print "Processing... folder", folderCount
            
            
    print "TOTAL FOLDERS CREATED:", folderCount
    



if __name__ == "__main__":
    
    
    print "STARTING!"
    
    #===============================================================================
    # Get input arguments, initialise class
    #===============================================================================
    parser = argparse.ArgumentParser(description= DESC_COMMENT) 

    # positional args
    parser.add_argument('rnaLibrary', metavar='rnaLibrary', type=str,
                         help='File path of catRAPID RNA library.')
    parser.add_argument('protLibrary', metavar='protLibrary', type=str,
                         help='File path of catRAPID Protein library.')
    parser.add_argument('rnaBatch', metavar='rnaBatch', type=int,
                         help='Number of RNAs to be batched on each process.')
    parser.add_argument('protBatch', metavar='protBatch', type=int,
                         help='Number of Proteins to be batched on each process.')
    parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
    parser.add_argument('catRAPIDTemplate', metavar='catRAPIDTemplate', type=str,
                         help='Template folder for with catRAPID library all vs all code. This should be the folder the script "run.library.both.sh".')    
    

    #gets the arguments
    args = parser.parse_args( ) 

    splitFile( args.rnaLibrary, args.rnaBatch, args.outputFolder)
    splitFile( args.protLibrary, args.protBatch, args.outputFolder)

    prepareFolders( args.rnaLibrary, args.protLibrary, args.outputFolder, args.catRAPIDTemplate)


    print "FINISHED!"
