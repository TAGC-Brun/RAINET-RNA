#!/usr/bin/python2.7
import os
import sys
import glob
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
# requires fastaq 3.11.0 to be installed

##### THIS SCRIPT SHOULD BE RUN FROM COMMAND LINE

###############################################################################
# output files will be written in same folder as input
# input fasta file must only contain one sequence
# inputFragmentsResult is output file from running catRAPID omics fragments (webserver automatically fragments if above 1200nt)
###############################################################################

inputFolder = sys.argv[1] #"/home/diogo/Documents/RAINET_data/catRAPID/webserver_results/HOTAIR_fragments/"
inputFasta = inputFolder + sys.argv[2] #"ENST00000424518_HOTAIR_001.fa"
inputFragmentsResult = inputFolder + sys.argv[3] #"ENST00000424518_HOTAIR_001.tsv"

###############################################################################

# function to read a output file from catRAPID omics (fragmenting RNA) and retrieve all the fragment coordinates
def read_input_fragments_result( input_fragments_result):
    # example
    #sp|Q9Y6V7|DDX49_HUMAN HOTAIR-001.cdna_1_253-404    -0.75    0.14    0.14    NO    yes    1    -    1.23
    #sp|A0AV96|RBM47_HUMAN HOTAIR-001.cdna_1_261-405    -0.12    0.50    0.99    NO    yes    1    -    1.61
    
    setOfCoordinates = set()
    
    with open( input_fragments_result, "r") as inFile:
        for line in inFile:           
            if line.startswith("#"):
                continue

            spl = line.split("\t")

            try:            
                spl2 = spl[0].split(" ")
                
                spl3 = spl2[1].split("_")
                
                coords = spl3[2]
                
                setOfCoordinates.add( coords)
                
            except IndexError as e:
                print e
                print "Failed to split line"
                print line
                print "expected format includes coodinates (fragmented transcript) such as: sp|A0AV96|RBM47_HUMAN HOTAIR-001.cdna_1_261-405    -0.12    0.50    0.99    NO    yes    1    -    1.61"
                raise Exception

    print "Found %s fragments" % len( setOfCoordinates)

    return setOfCoordinates


# Launches commands to clip fasta file based on given set of coordinates
def clip_fasta(input_fasta, set_of_coordinates):

    #e.g. fastaq add_indels --delete HOTAIR:1-219 --delete HOTAIR:416-9999999 test.fa tgest.fa

    ### replace ":" in fasta file to "_"
    modInputFile = input_fasta.replace(".fa","") + ".modified"
    cmd = "sed 's/:/_/g' %s | sed 's/ /_/g' > %s" % ( input_fasta, modInputFile)
#    os.system(cmd)
    SubprocessUtil.run_command( cmd)

    ### Get sequence header
    name = ""
    with open( modInputFile, "r") as inFile:
        for line in inFile:
            if ">" in line:
                name = line.strip()[1:]
                break

    if len( name) < 1:
        print "Failure to detect name of fasta sequence"
        raise Exception


    ### Clip fasta sequence for each set of coordinates
    for coord in set_of_coordinates:
                
        start, end = coord.split("-")

        outFile = input_fasta + "_" + start + "-" + end

        # delete for start of sequence to start of wanted coordinate
        tag1 = name + ":1-" + start
        # deleted from end of wanted coordinate to end of sequence
        tag2 = name + ":" + end + "-999999"

        command = "fastaq add_indels --delete %s --delete %s %s %s" % ( tag1, tag2, modInputFile, outFile)

        SubprocessUtil.run_command( command)
#        os.system( command)

    ### Change fasta headers 

    #e.g. fastaq enumerate_names ENST00000424518_HOTAIR_001.fa_2288-2422 ENST00000424518_HOTAIR_001.fa_2288-2422_renamed --suffix ENST00000424518_HOTAIR_001.fa_2288-2422    

    newFastas = glob.glob( inputFolder + "/*.fa_*")
    
    for fasta in newFastas:
        fileName = fasta.split("/")[-1]
        command = "fastaq enumerate_names %s %s --suffix %s" % ( fileName, fileName + "_renamed", fileName)
        SubprocessUtil.run_command( command)
        
    ### Concatenate files into a multi fasta
    newFileName = input_fasta + ".all_fragments.fa"
    command = "cat *renamed* > %s" % ( newFileName)    
    SubprocessUtil.run_command( command)

#     ### convert into oneline format (for catRAPID)
#     command = "fastaq to_fasta -l 0 %s %s" % ( newFileName, newFileName + ".nospace")
#     SubprocessUtil.run_command( command)
#     command = "sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' %s | sed 's/>/\n>/g' > %s" % ( newFileName + ".nospace", newFileName + ".oneline")
#     SubprocessUtil.run_command( command)

    ### clean up
    command = "rm *.fa_*"
    SubprocessUtil.run_command( command)



setOfCoordinates = read_input_fragments_result( inputFragmentsResult)

clip_fasta( inputFasta, setOfCoordinates)

print "FINISHED!"