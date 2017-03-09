
import sys
import os

# 09-Mar-2017 Diogo Ribeiro
# Routine to merge and process RPISeq raw results into a format readable to my scripts (catRAPID format)

inputFiles = sys.argv[1] # file with list of RPISeq output files to parse here, full paths must be given
outputFolder = sys.argv[2] # folder where to write output file
proteinID = sys.argv[3] # ID of the interacting protein, to appear in output file

### Input example (RPISeq output)
# RNA ID    RF Classifier    SVM Classifier
# >ENST00000437898    0.75    0.541
# >ENST00000452728    0.75    0.774
# >ENST00000422049    0.8    0.451
# (split into several files, each one with a header)
# the protein name is defined beforehand, e.g. P52298
# RPISeq provides two scores, one output file will be created for each

### Wanted Output
# sp|P52298|NCBP2_HUMAN ENST00000625358   40.15   0.88    0.08
# sp|P52298|NCBP2_HUMAN ENST00000625632   49.31   0.94    0.29
# sp|P52298|NCBP2_HUMAN ENST00000627767   48.38   0.94    0.29
# (no header)

outFileRF = open( outputFolder + "/RIPSeq_RF_formatted.out", "w")
outFileSVM = open( outputFolder + "/RIPSeq_SVM_formatted.out", "w")

with open( inputFiles, "r") as initialFile:
    for li in initialFile:
        # open each file
        with open( li.strip(), "r") as inFile:
            header = inFile.readline()
            for line in inFile:
                spl = line.strip().split()
                
                if spl < 3:
                    print "Problem with input file: %s, %s" % (inFile, line)

                txID = spl[0].replace(">","")
                rfScore = spl[1]
                svmScore = spl[2]

                # write RF scores                
                outFileRF.write( "sp|%s|sp %s\t%s\n" % (proteinID, txID, rfScore))
                
                # write SVM scores
                outFileSVM.write( "sp|%s|sp %s\t%s\n" % (proteinID, txID, svmScore))

outFileRF.close()
outFileSVM.close()



            