

from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
import subprocess
import os
import glob

#===============================================================================
# 12-Fev-2016 Diogo Ribeiro
# Script to process GTEx tissue expression files
#
# Objective is to produce basic statistics on the expression levels of transcripts,
# normalise and modify the data for Rainet database insertion. I.e. reduce sample dimension
# E.g. have a single expression value per transcript per tissue.
#===============================================================================

#===============================================================================
# Plan
# 
# Read file GTEx_Data_V6_Annotations_SampleAttributesDS.txt, make dictionary of sample IDs -> tissue
# Make basic stats on number of samples per tissue
# Read header of GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm, make column indexes for each tissue (e.g. in the 8k samples, which ones belong to blood, which to heart etc)
# Read GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm entries and group results of samples of same tissue.
# Measure variability of expression per tissue, and variability of expression per RNA, look at GTEx papers as well
# Based on distribution of values, merge data so that we have a single RPKM value for each transcript-tissue pair
# Add data to RAINET database!
#===============================================================================

class ProcessGTExData( object ):
    
    TISSUE_ANNOTATIONS = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/RNA/expression_test/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
    TISSUE_EXPRESSION = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/RNA/expression_test/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm_test.txt"

    TISSUE_ANNOTATIONS_KEY = "SAMPID"
    TISSUE_ANNOTATIONS_VALUE = "SMTS" # Which level of tissue definition to use for grouping samples, could also be SMTSD

    # Python files
    ANNOTATION_OUTPUT_FILE = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/dataFiles/annotation.csv"
    EXPRESSION_OUTPUT_FILE = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/dataFiles/expression.csv"

    # R files    
    WORKING_DIR = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/expressionStatistics"
    SWEAVE_R_SCRIPT = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/expressionStatistics/GTEx_statistics.R"

    # Files to be used by R
    R_REQUIRED_FILES = [SWEAVE_R_SCRIPT, WORKING_DIR, ANNOTATION_OUTPUT_FILE]
    
    # Files to be removed before each run of this script
    FILES_TO_REFRESH = [ANNOTATION_OUTPUT_FILE, EXPRESSION_OUTPUT_FILE, WORKING_DIR+"/*.tex", WORKING_DIR+"/*.log", WORKING_DIR+"/*.aux", WORKING_DIR+"/*.pdf"]


    def read_annotations(self):
        
        inHandler = FileUtils.open_text_r( ProcessGTExData.TISSUE_ANNOTATIONS)
        
        # Read header
        headerLine = inHandler.readline()
        header = headerLine.split("\t")
            
        # Build the map of parameter to header index
        # If a parameter name is not in the header raise exception
        headerIndexMap = {}
        for parameter in header:
            try:
                index = header.index( parameter)
            except ValueError:
                raise RainetException( "read_annotations : The parameter '" + parameter + "' is not in the header: ")

            headerIndexMap[ parameter] = index


        # Read the rest of the file
        counter = 0
        sampleTissue = {} # key -> sample ID, val -> tissue
        tissueSample = {} # key -> tissue, val -> sample ID
        for line in inHandler:
            if line != '' and not line.startswith( "#") and len(line) > 0:
                spl = line.split( "\t")

                #Check if the line has as many columns as header
                if len( spl) > len( header):
                    raise RainetException( "read_annotations : Bad line format in file. More columns than in header: " + str(len(spl)) )

                keyIdx = headerIndexMap[ProcessGTExData.TISSUE_ANNOTATIONS_KEY]
                valIdx = headerIndexMap[ProcessGTExData.TISSUE_ANNOTATIONS_VALUE]

                sample = spl[keyIdx]
                tissue = spl[valIdx]

                if len( tissue) == 0 or tissue == None:
                    Logger.get_instance().warning( "read_annotations : Sample without tissue information: "+ sample + ". Skipping this sample.." )
                    continue

                if sample not in sampleTissue:
                    sampleTissue[ sample] = tissue
                else:
                    raise RainetException( "read_annotations : Duplicate sargumentFilesample entry: "+ sample )
                    
                if tissue not in tissueSample:
                    tissueSample[ tissue] = []
                tissueSample[ tissue].append( sample)

                counter+= 1

        if counter != len( sampleTissue):
            Logger.get_instance().warning( "read_annotations : Number of lines read different from total number of samples: " + str( counter) + " vs " + str( len( sampleTissue)) )

        Logger.get_instance().info( "read_annotations : Number of samples found = "+ str(len( sampleTissue) ) )
        Logger.get_instance().info( "read_annotations : Number of tissues found = "+ str(len( tissueSample) ) )

        inHandler.close()

        outHandler = FileUtils.open_text_w( ProcessGTExData.ANNOTATION_OUTPUT_FILE )
        
        for tissue in tissueSample:
            outHandler.write( "%s,%s\n" % ( tissue, len(tissueSample[tissue])) )
        outHandler.close()


    def run_statistics(self):

        # confirm that required input files are present
        for filePath in ProcessGTExData.R_REQUIRED_FILES:
            if not os.path.exists(filePath):
                raise RainetException( "run_statistics : Input file is not present: " + filePath )
                
        # launch the analysis
        command = "cd " + os.path.dirname(ProcessGTExData.SWEAVE_R_SCRIPT) + "; Rscript %s %s %s" % ( ProcessGTExData.SWEAVE_R_SCRIPT, ProcessGTExData.WORKING_DIR, ProcessGTExData.ANNOTATION_OUTPUT_FILE)
                
        Logger.get_instance().info( "run_statistics : Running command : "+command)
        
        # run the command in a subprocess
        outfile = open( ProcessGTExData.WORKING_DIR+"/analysis.log", "wb")
        p = subprocess.Popen( command, shell=True, stdout=outfile, stderr=outfile)
        while True:
            if p.poll() != None:
                break


    # #
    # Method to remove files that are recreated each time.
    def clean_files(self):
        
        for filePath in ProcessGTExData.FILES_TO_REFRESH:
            files = glob.glob(filePath)
            for fi in files:
                os.remove(fi)


if __name__ == "__main__":
    
    try:

        # Create Logger instance by using the first log action.
        Logger.get_instance().info( "ProcessGTExData : Starting..." )

        # Start chrono
        Timer.get_instance().start_chrono()

        run = ProcessGTExData()

        run.clean_files()

        Timer.get_instance().step( "reading annotation file..")    
        run.read_annotations()
        Timer.get_instance().step( "reading expression file..")    
 
 
        run.run_statistics()


    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of ProcessGTExData. Aborting :\n" + rainet.to_string())

    # Stop the chrono      
    Timer.get_instance().stop_chrono( "ProcessGTExData : Finished" )

