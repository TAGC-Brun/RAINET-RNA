
import subprocess
import os
import glob
import numpy as np
from scipy import stats

from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer

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

    # GTEx Input files
    TISSUE_ANNOTATIONS = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/RNA/expression_test/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
    TISSUE_EXPRESSION = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/RNA/expression_test/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm_test.txt"

    # Fields from TISSUE_ANNOTATIONS file that we want to keep
    TISSUE_ANNOTATIONS_KEY = "SAMPID"
    TISSUE_ANNOTATIONS_VALUE = "SMTS" # Which level of tissue definition to use for grouping samples, could also be SMTSD

    # Fields from TISSUE_EXPRESSION that are special
    TISSUE_EXPRESSION_SPECIAL_COLUMNS = [ "TargetID", "Gene_Symbol", "Chr", "Coord"]

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


    def __init__(self):
        pass
        

    # #
    # Read GTEx annotations file and build mapping between tissue and samples
    def read_tissue_annotations(self):
        
        inHandler = FileUtils.open_text_r( ProcessGTExData.TISSUE_ANNOTATIONS)
        
        # Read header
        headerLine = inHandler.readline()
        header = headerLine.split("\t")
            
        # Build the map of parameter to header index
        headerIndexMap = {}
        for parameter in header:
            index = header.index( parameter)
            headerIndexMap[ parameter] = index

        # Read the rest of the file
        counter = 0
        sampleTissue = {} # key -> sample ID, val -> tissue
        tissueSample = {} # key -> tissue, val -> sample ID
        problematicSamples = set() # stores sample IDs that are problematic/erronic
        for line in inHandler:
            if line != '' and not line.startswith( "#") and len(line) > 0:
                spl = line.split( "\t")

                #Check if the line has as many columns as header
                if len( spl) > len( header):
                    raise RainetException( "read_tissue_annotations : Bad line format in file. More columns than in header: " + str(len(spl)) )

                keyIdx = headerIndexMap[ProcessGTExData.TISSUE_ANNOTATIONS_KEY]
                valIdx = headerIndexMap[ProcessGTExData.TISSUE_ANNOTATIONS_VALUE]

                sample = spl[keyIdx]
                tissue = spl[valIdx]

                if len( tissue) == 0 or tissue == None:
                    Logger.get_instance().warning( "read_tissue_annotations : Sample without tissue information: "+ sample + ". Skipping this sample.." )
                    problematicSamples.add( sample)
                    continue

                if sample not in sampleTissue:
                    sampleTissue[ sample] = tissue
                else:
                    raise RainetException( "read_tissue_annotations : Duplicate sargumentFilesample entry: "+ sample )
                    
                if tissue not in tissueSample:
                    tissueSample[ tissue] = []
                tissueSample[ tissue].append( sample)

                counter+= 1

        if counter != len( sampleTissue):
            Logger.get_instance().warning( "read_tissue_annotations : Number of lines read different from total number of samples: " + str( counter) + " vs " + str( len( sampleTissue)) )

        Logger.get_instance().info( "read_tissue_annotations : Number of samples found = "+ str(len( sampleTissue) ) )
        Logger.get_instance().info( "read_tissue_annotations : Number of tissues found = "+ str(len( tissueSample) ) )

        inHandler.close()

        # Write file to be read by R
        outHandler = FileUtils.open_text_w( ProcessGTExData.ANNOTATION_OUTPUT_FILE )       
        #outHandler.write("%s,%s\n" % (ProcessGTExData.TISSUE_ANNOTATIONS_VALUE,"'#_samples'") )
        for tissue in tissueSample:
            outHandler.write( "%s,%s\n" % ( tissue, len(tissueSample[tissue])) )
        outHandler.close()

        return sampleTissue, tissueSample, problematicSamples

    # #
    # 
    def read_transcript_expression(self, sampleTissue, tissueSample, problematicSamples):

        inHandler = FileUtils.open_text_r( ProcessGTExData.TISSUE_EXPRESSION)
        
        # Read header
        headerLine = inHandler.readline()
        header = headerLine.split("\t")
        lenHeader = len( header)

        # Build the map of parameter to header index
        # If a sample name is not in the header raise exception
        headerIndexMap = {} # key -> tissue , val -> list of indexes
        processedSamples = 0
        for i in xrange( lenHeader ):
            sample = header[i].strip()
            if sample not in ProcessGTExData.TISSUE_EXPRESSION_SPECIAL_COLUMNS:
                if sample not in sampleTissue:
                    if sample not in problematicSamples:
                        raise RainetException( "read_transcript_expression : sample ID not found in tissue annotations: "+ sample )
                    else:
                        Logger.get_instance().warning( "read_transcript_expression : Known sample without tissue info present in header: "+ sample + ". Will ignore information for this sample." )
                else:
                    tiss = sampleTissue[sample]
                    if tiss not in headerIndexMap:
                        headerIndexMap[tiss] = set()

                    headerIndexMap[tiss].add( i)
                    processedSamples+= 1


        summ = 0
        for tiss in headerIndexMap:
            summ+= len(headerIndexMap[tiss])
#             # trying to use ranges instead of actual coordinates. Not possible because of erroneous data
#             minIndex = min(headerIndexMap[ tiss])
#             maxIndex = max(headerIndexMap[ tiss])
#             for i in xrange( minIndex, maxIndex):
#                 if i not in headerIndexMap[ tiss]:
#                     print ( tiss, i, minIndex, maxIndex)  
        if summ != processedSamples:
            raise RainetException( "read_transcript_expression : index insertion was not performed correctly: ")

        # initialise dictionaries which will contain results
        txExpressionTissue = {} #transcript expression per tissue. key -> transcript ID, val -> dict: key -> tissue, val -> list of expression values
        expressionTissue = {} #expression per tissue. key -> tissue, val -> list of expression values (i.e. across transcripts)
        for tiss in headerIndexMap:
            expressionTissue[tiss] = []

        # Read each non-header line of expression file
        count = 0
        for line in inHandler:

            spl = line.strip().split("\t")
            
            if len( spl) != lenHeader:
                raise RainetException( "read_transcript_expression : number of items in line is not the expected: "+ str( len( spl)) )

            # Getting transcript ID
            #
            # Several hundred transcript IDs are the merge of two transcript IDs. 
            # Testing a few, there was always one that was deprecated and so that seems to be the reason for this.
            # I decided to include both transcripts in analysis, both sharing the same expression values            
            tx = spl[0]
            txSpl = tx.split("_")
            if len( txSpl) > 2: 
                raise RainetException( "read_transcript_expression : transcript ID is not expected: "+ tx )

            for tx in txSpl:
                # GTEx should use only Ensembl transcript IDs (ENST*), raise exception if something else is found
                if not tx.startswith("ENST"):
                    raise RainetException( "read_transcript_expression : transcript ID is not expected: "+ tx )
    
                # GTEx transcript IDs have a different termination which is not used in Ensembl (e.g. "ENST00000002501.6"). This is being removed here.
                if "." not in tx:
                    transcriptID = tx
                else:
                    txSpl = tx.split(".")
                    if len( txSpl) == 2:
                        transcriptID = txSpl[0]
                    else:
                        raise RainetException( "read_transcript_expression : transcript ID is not expected: "+ tx )

                # initialise dict for transcript
                if transcriptID not in txExpressionTissue:
                    txExpressionTissue[transcriptID] = {}
    
                for tiss in headerIndexMap:
    
                    # initialise for transcript-tissue pair
                    if tiss not in txExpressionTissue[transcriptID]:
                        txExpressionTissue[transcriptID][tiss] = []
                    
                    # add data to data structures
                    for idx in headerIndexMap[tiss]:
                        try:
                            value = float(spl[idx])
                        except ValueError:
                            raise RainetException( "read_transcript_expression : expression value is non-numeric: "+ spl[idx] )                           
                        txExpressionTissue[transcriptID][tiss].append( value)
                        expressionTissue[tiss].append( value)

            count+= 1
            if count % 1000 == 0:
                Logger.get_instance().info( "read_transcript_expression : reading file.. %s lines done." % count)


        # #
        # Write into files for R

        outHandler = FileUtils.open_text_w( ProcessGTExData.EXPRESSION_OUTPUT_FILE )        

        # File with data per tissue
        for tiss in expressionTissue:
            line = tiss+","
            for val in expressionTissue[tiss]:
                line+= str(val)+","
            line = line[:-1]+"\n" # remove last comma
            outHandler.write(line)
            
        outHandler.close()

        print ("Total tissues", len(expressionTissue))
        print ("Total transcripts", len(txExpressionTissue)) #check that total transcripts is same as gencode v19 # I should have it as unittest


        ###make same report as before but with the samples that are on expression file

        #e.g. validation
#         grep ENST00000002501.6 GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm_test.txt | cut -f7927
#         1.721180


    # #
    # Run Rscript to produce Sweave file and consequent pdf report, using the data written by this script
    def run_statistics(self):

        # confirm that required input files are present
        for filePath in ProcessGTExData.R_REQUIRED_FILES:
            if not os.path.exists(filePath):
                raise RainetException( "run_statistics : Input file is not present: " + filePath )
                
        # launch the analysis
        command = "cd " + os.path.dirname(ProcessGTExData.SWEAVE_R_SCRIPT) + "; Rscript %s %s %s %s" % ( ProcessGTExData.SWEAVE_R_SCRIPT, ProcessGTExData.WORKING_DIR, ProcessGTExData.ANNOTATION_OUTPUT_FILE, ProcessGTExData.EXPRESSION_OUTPUT_FILE)

        Logger.get_instance().info( "run_statistics : Running command : "+command)
        
        # run the command in a subprocess
        p = subprocess.Popen( command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while True:
            if p.poll() != None:
                break

        stderrText = p.stderr.read().decode("UTF-8").strip()
        stdoutText = p.stdout.read().decode("UTF-8").strip()
    
        if len(stderrText) > 0:
            Logger.get_instance().warning("run_statistics : STDERR:\n"+stderrText)

        if p.returncode != 0:
            Logger.get_instance().error("run_statistics : ERROR: Return code:"+str(p.returncode)+"\t"+command)
        
        if len(stdoutText) > 0:
            Logger.get_instance().info("run_statistics : STDOUT:\n"+stdoutText)



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

        # Initialise class
        run = ProcessGTExData()

        # Wipe output and log from previous runs
        run.clean_files()

        # Read annotation file
        Timer.get_instance().step( "reading annotation file..")    
        sampleTissue, tissueSample, problematicSamples = run.read_tissue_annotations()
        
        # Read expression file, using annotations
        Timer.get_instance().step( "reading expression file..")    
        run.read_transcript_expression(sampleTissue, tissueSample, problematicSamples)
        
        # Run R scripts / report
        run.run_statistics()

    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of ProcessGTExData. Aborting :\n" + rainet.to_string())

    # Stop the chrono      
    Timer.get_instance().stop_chrono( "ProcessGTExData : Finished" )

