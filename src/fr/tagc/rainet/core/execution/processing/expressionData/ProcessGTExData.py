
import sys
import os
import argparse
import glob
import numpy as np
import random 

from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil

#===============================================================================
# Started 12-Fev-2016 
# Diogo Ribeiro
# Script to process GTEx tissue expression files
#
# Objective is to produce basic statistics on the expression levels of transcripts,
# normalise and modify the data (e.g. reduce sample dimension) for Rainet database insertion. 
#===============================================================================

#===============================================================================
# General plan:
# 
# - Read file GTEx_Data_V6_Annotations_SampleAttributesDS.txt, make dictionary of sample IDs -> tissue
# - Read header of GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm, make column indexes for each tissue (e.g. in the 8k samples, which ones belong to blood, which to heart etc)
# - Read GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm entries and group results of samples of same tissue.
# - Produce report on variability of expression per tissue, and variability of expression per RNA
# - Average data so that we have a single RPKM value for each transcript-tissue pair
# - Write output file to be used RAINET database insertion
#===============================================================================

#===============================================================================
# Processing notes:
# 
# 0 - This module has a unittest test which should be ran before this script
# 1 - Some sample entries in GTEx_Data_V6_Annotations_SampleAttributesDS.txt
#     do not contain tissue ("SMTSD") information, these are excluded
# 2 - Some expression entries (~700 out of 200k~) in GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt
#     have a merged transcriptID identifier (e.g. ENSTXXXXXX_ENSTXXXXXX), whereas one
#     of the identifiers points to a deprecated transcriptID in a later Ensembl version (e.g. v83).
#     We keep both identifiers and attribute the same expression values to each.
#
#===============================================================================


class ProcessGTExData( object ):

    # Fields from TISSUE_ANNOTATIONS file that we want to keep
    TISSUE_ANNOTATIONS_KEY = "SAMPID"

    # Fields from TISSUE_EXPRESSION that are not samples
    TISSUE_EXPRESSION_SPECIAL_COLUMNS = [ "TargetID", "Gene_Symbol", "Chr", "Coord"]

    # Python output report files
    ANNOTATION_OUTPUT_FILE = "annotation.csv"
    EXPRESSION_OUTPUT_FILE = "expression.csv"
    TX_EXPRESSION_OUTPUT_FOLDER = "transcript_expression/"
    TX_EXPRESSION_AVG_OUTPUT_FILE = "transcript_expression_metrics.tsv"
    TX_EXPRESSION_AVG_OUTPUT_FILE_NO_OUTLIERS = "transcript_expression_metrics_no_outliers.tsv"
    
    # R script files    
    WORKING_DIR = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/expressionStatistics"
    SWEAVE_R_SCRIPT = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/expressionData/expressionStatistics/GTEx_statistics.R"

    # Files to be used by R
    R_REQUIRED_FILES = [ANNOTATION_OUTPUT_FILE, EXPRESSION_OUTPUT_FILE, TX_EXPRESSION_AVG_OUTPUT_FILE, TX_EXPRESSION_AVG_OUTPUT_FILE_NO_OUTLIERS]
    
    # Predefined tissues
    PREDEFINED_TISSUES = ["Muscle - Skeletal","Whole Blood","Stomach","Minor Salivary Gland","Kidney - Cortex"]

    # Sampling for report
    RNA_SAMPLE_NUMBER = 50 # for the tissue expression boxplots #note that final picked number may be slightly different
    RNA_TEMP_FILE = WORKING_DIR+"/RNA_TEMP_FILE" # temporary file for the individual tx plots
    RNA_SAMPLE_NUMBER_PER_TX = 12 # number of tx used for the individual tx report


    def __init__(self, annotation_file, expression_file, tissue_level, minimum_samples, output_folder, write_report ):

        # assign attributes from system arguments
        self.annotationFile = annotation_file
        self.expressionFile = expression_file
        self.tissueLevel = tissue_level
        self.minimumSamples = minimum_samples
        self.outputFolder = output_folder
        self.writeReport = write_report

    # #
    # Read GTEx sample annotations file and build mapping between tissue/body-parts and sample IDs
    #
    # @return sampleTissue dictionary containing correspondence of sample IDs to tissue name
    # @return tissueSample dictionary containing correspondence of tissue to list of sample IDs
    # @return problematicSamples contains set of sample IDs that contained erroneous data and will be excluded
    def read_tissue_annotations(self):
        
        inHandler = FileUtils.open_text_r( self.annotationFile)
        
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
        tissueSample = {} # key -> tissue, val -> list of sample IDs
        problematicSamples = set() # stores sample IDs that are problematic/erronic
        for line in inHandler:
            if line != '' and not line.startswith( "#") and len( line) > 0:
                spl = line.split( "\t")

                #Check if the line has as many columns as header
                if len( spl) > len( header):
                    raise RainetException( "read_tissue_annotations : Bad line format in file. More columns than in header: " + str(len(spl)) )

                keyIdx = headerIndexMap[ ProcessGTExData.TISSUE_ANNOTATIONS_KEY]
                valIdx = headerIndexMap[ self.tissueLevel]

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

        return sampleTissue, tissueSample, problematicSamples


    # #
    # Read GTEx RPKM expression file 
    # uses previous sample-tissue mapping to group data on transcript-tissue
    #
    # @param sample_tissue : dictionary - contains correspondence of sample IDs to tissue name
    # @param tissue_sample : dictionary - contains correspondence of tissue to list of sample IDs
    # @param problematic_samples : set - contains sample IDs that contained erroneous data and will be excluded
    # 
    # Writes several files for GTEx pdf report
    #
    # @return txExpressionTissue : dictionary - contains expression data for each transcript, over all tissues and all samples
    def read_transcript_expression(self, sample_tissue, tissue_sample, problematic_samples):

        inHandler = FileUtils.open_text_r( self.expressionFile)
        
        #===============================================================================
        # Read header, map sample ID to tissue names, remove tissues with less than minimumSamples
        #===============================================================================
        headerLine = inHandler.readline()
        header = headerLine.split( "\t")
        lenHeader = len( header)

        # Build the map of parameter to header index
        # If a sample name is not in the header raise exception
        headerIndexMap = {} # key -> tissue , val -> list of indexes
        processedSamples = 0
        for i in xrange( lenHeader ):
            sample = header[i].strip()
            if sample not in ProcessGTExData.TISSUE_EXPRESSION_SPECIAL_COLUMNS:
                if sample not in sample_tissue:
                    if sample not in problematic_samples:
                        raise RainetException( "read_transcript_expression : sample ID not found in tissue annotations: "+ sample )
                    else:
                        Logger.get_instance().warning( "read_transcript_expression : Known sample without tissue info present in header: "+ sample + ". Will ignore information for this sample." )
                else:
                    tiss = sample_tissue[ sample]
                    if tiss not in headerIndexMap:
                        headerIndexMap[ tiss] = set()

                    headerIndexMap[ tiss].add( i)
                    processedSamples+= 1

        # Exclude tissues with less than provided minimum number of samples
        summ = 0
        toDelete = set()
        for tiss in headerIndexMap:
            lenSamples = len( headerIndexMap[ tiss])            
            if lenSamples < self.minimumSamples:
                toDelete.add( tiss)

            summ+= lenSamples
            
        if summ != processedSamples:
            raise RainetException( "read_transcript_expression : index insertion was not performed correctly: ")
        
        for tiss in toDelete:
            Logger.get_instance().info( "read_transcript_expression : tissue removed from minimumSamples filter: " + tiss)
            del headerIndexMap[ tiss]


        #===============================================================================
        # Read the rest of the file, i.e. the RPKM values
        #===============================================================================
        
        # initialise dictionaries which will contain results
        txExpressionTissue = {} #transcript expression per tissue. key -> transcript ID, val -> dict: key -> tissue, val -> list of expression values
        expressionTissue = {} #expression per tissue. key -> tissue, val -> list of expression values (i.e. across transcripts)
        for tiss in headerIndexMap:
            expressionTissue[ tiss] = []

        # pick random index of transcripts which will be used for report
        # get number of lines in file, quickest way possible
        wcResult = SubprocessUtil.run_command("wc -l %s" % ( self.expressionFile), return_stdout = 1, verbose = 0 )
        lengthOfFile = int( wcResult.split(" ")[0])
        txSampling = random.sample( xrange( 0, lengthOfFile), ProcessGTExData.RNA_SAMPLE_NUMBER)
        
        # loop per line
        count = 0
        countValues = 0
        for line in inHandler:

            spl = line.strip().split( "\t")
            
            if len( spl) != lenHeader:
                raise RainetException( "read_transcript_expression : number of items in line is not the expected: "+ str( len( spl)) )

            # Getting transcript ID
            #
            # Several hundred transcript IDs are the merge of two transcript IDs. 
            # Testing a few, there was always one that was deprecated and so that seems to be the reason for this.
            # I decided to include both transcripts in analysis, both sharing the same expression values            
            tx = spl[0]
            txSpl = tx.split("_")
 
            for tx in txSpl:
                # GTEx should use only Ensembl transcript IDs (ENST*), raise exception if something else is found
                if not tx.startswith( "ENST"):
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
                    txExpressionTissue[ transcriptID] = {}
     
                # loop over each tissue to retrieve data for all sample indexes of that tissue
                for tiss in headerIndexMap:
     
                    # add data to local data structures
                    # keep data in temporary list to avoid excessive memory requirements, then turn into numpy array
                    localList = [] 
                    for idx in headerIndexMap[ tiss]:
                        try:
                            value = float( spl[ idx])
                        except ValueError:
                            raise RainetException( "read_transcript_expression : expression value is non-numeric: "+ spl[ idx] )                           
 
                        localList.append( value)
                        countValues+= 1
 
                        # Whether this is running in the test file or not, get only values for a limited amount of transcripts             
                        if count in txSampling:
                            expressionTissue[ tiss].append( value)
     
                    # initialise transcript-tissue pair
                    if tiss not in txExpressionTissue[ transcriptID]:
                        txExpressionTissue[ transcriptID][ tiss] = np.array( localList)


            count+= 1
            if count % 1000 == 0:
                Logger.get_instance().info( "read_transcript_expression : reading file.. %.2f%% done." % ( count * 100.0 / lengthOfFile ) )
#             if count > 10000: #to remove
#                 break

        Logger.get_instance().info( "read_transcript_expression : Total tissues processed = " + str( len( headerIndexMap)) )
        Logger.get_instance().info( "read_transcript_expression : Total transcripts processed = " + str( len( txExpressionTissue)) )
        Logger.get_instance().info( "read_transcript_expression : Total expression values processed = " + str( countValues ) )
        
        #===============================================================================
        # Write variables into files for report
        #===============================================================================
        
        # File with number of samples per tissue
        #
        outHandler = FileUtils.open_text_w( self.outputFolder + ProcessGTExData.ANNOTATION_OUTPUT_FILE )       
        for tissue in headerIndexMap:
            outHandler.write( "%s,%s\n" % ( tissue, len(headerIndexMap[tissue])) )
        outHandler.close()

        
        # Files with data per tissue for sample of transcripts: Aim is to use this data to produce a boxplot for each tissue
        #
        # E.g. of wanted output: tissue and all its values on each line
        # Thyroid,99.987946,51.082527,73.664589,107.893951
        # Kidney,78.215836,94.459152,104.612236,143.620743,91.351257,134.171814,144.965485,111.592468
        # Pancreas,92.215179,98.487106

        outHandler = FileUtils.open_text_w( self.outputFolder + ProcessGTExData.EXPRESSION_OUTPUT_FILE )        
        
        for tiss in expressionTissue:
            line = tiss+","
            for val in expressionTissue[tiss]:
                line+= str(val)+","
            line = line[:-1]+"\n" # remove last comma
            outHandler.write(line)
            
        outHandler.close()


        # File with data per transcript: aim is to produce a histogram of expression values of a few randomly selected transcripts
        #
        # Create folder to contain several files to be read by R, one file per transcript.

        if not os.path.exists( self.outputFolder + ProcessGTExData.TX_EXPRESSION_OUTPUT_FOLDER):
            os.mkdir( self.outputFolder + ProcessGTExData.TX_EXPRESSION_OUTPUT_FOLDER)

        txSample = random.sample( txExpressionTissue.keys(), ProcessGTExData.RNA_SAMPLE_NUMBER_PER_TX)
        
        for tx in txSample:
            outHandler = FileUtils.open_text_w( self.outputFolder + ProcessGTExData.TX_EXPRESSION_OUTPUT_FOLDER+"/"+tx )
            txData = txExpressionTissue[tx]
            for tiss in txData:
                if tiss in ProcessGTExData.PREDEFINED_TISSUES:
                    line = tiss+","
                    for val in txData[tiss]:
                        line+= str(val)+","
                    line = line[:-1]+"\n" # remove last comma
                    outHandler.write(line)

            outHandler.close()


        return txExpressionTissue


    # #
    # Function to produce single value for each tissue-transcript pair
    #
    # Write to final file for insertion into RAINET database
    #
    #@param tx_expression_tissue : dictionary produced by read_transcription_expression
    def average_sample_values(self, tx_expression_tissue):
        
        # Outputs TSV file
        # e.g. TranscriptID    TissueName      ExprMean        ExprStd ExprMedian      CoefVariation   Max
        #      ENST00000496116 Muscle - Skeletal       0.101   0.103   0.069   1.02415982815   0.426

        # open output files        
        outHandler = FileUtils.open_text_w( self.outputFolder + ProcessGTExData.TX_EXPRESSION_AVG_OUTPUT_FILE )
        outHandlerNoOutliers = FileUtils.open_text_w( self.outputFolder + ProcessGTExData.TX_EXPRESSION_AVG_OUTPUT_FILE_NO_OUTLIERS )

        # Header
        outHandler.write( "transcript_id\ttissue_name\trpkm_mean\trpkm_std\trpkm_median\tcoef_variation\tmax\n")
        outHandlerNoOutliers.write( "transcript_id\ttissue_name\trpkm_mean\trpkm_std\trpkm_median\tcoef_variation\tmax\n")
        
        percOfRemoved = []
        
        # Calculate several metrics of RPKM expression
        count = 0
        for tx in tx_expression_tissue:
            count+=1
            if count % 10000 == 0:
                Logger.get_instance().info( "read_transcript_expression : writing to file.. %s lines done." % count)

            for tiss in tx_expression_tissue[ tx]:

                # get array of values for transcript-tissue
                sampleArray = tx_expression_tissue[ tx][ tiss]

                # #
                # Write file without removing outliers

                # compute descriptive statistics
                maxi = np.max( sampleArray)
                std = np.std( sampleArray) 
                mean = np.mean( sampleArray)
                median = np.median( sampleArray)

                if np.isnan(std) or np.isnan(mean) or mean == 0:
                    coefVar = 0
                else:
                    coefVar = float(std) / float(mean) 

                outHandler.write("%s\t%s\t%.3f\t%.3f\t%.3f\t%s\t%.3f\n" % ( tx, tiss, mean, std, median, coefVar, maxi ) )

                # #
                # Write file removing outliers
                
                valuesInRange = self.remove_outliers( sampleArray)

                numberRemoved =  len(sampleArray) - len(valuesInRange)
                percRemoved = numberRemoved * 100.0 / len(sampleArray)                
                percOfRemoved.append(percRemoved)
        
                # print ( "remove_outliers : Number of samples removed: %s out of %s (%.2f%%) " % \
                #          ( numberRemoved, len(sampleArray), percRemoved ) )
                
                # compute descriptive statistics
                maxi = np.max( valuesInRange)
                std = np.std( valuesInRange) 
                mean = np.mean( valuesInRange)
                median = np.median( valuesInRange)

                if np.isnan(std) or np.isnan(mean) or mean == 0:
                    coefVar = 0
                else:
                    coefVar = float(std) / float(mean) 

                outHandlerNoOutliers.write("%s\t%s\t%.3f\t%.3f\t%.3f\t%s\t%.3f\n" % ( tx, tiss, mean, std, median, coefVar, maxi ) )

        Logger.get_instance().info("read_transcript_expression : Mean and median %% of samples removed with outlier removal:\t%.2f%%\t%.2f%%" % (np.mean(percOfRemoved),np.median(percOfRemoved)) )

        # remove large dictionary from memory after writing the files
        del tx_expression_tissue

    # #
    # Method to exclude outlier values from a numpy array
    # Remove outliers by excluding values out of the range: Q1 - 1.5xIQR : Q3 + 1.5xIQR
    #
    # @param : numpy array - list of values
    # @return : numpy array - list of values without outliers
    def remove_outliers(self, array):

        # Note: Q1 = 25th percentile, Q3 = 75th percentile
        q1, q3 = np.percentile( array, [ 25, 75])
        iqr = q3 - q1
        rangeMin = q1 - 1.5 * iqr
        rangeMax = q3 + 1.5 * iqr 
        valuesInRange = np.array( [x for x in array if rangeMin <= x <= rangeMax])

        return valuesInRange


    # #
    # Run Rscript to produce Sweave file and consequent pdf report, using the data written by this script
    def run_statistics(self):

        # confirm that required input files are present
        for filePath in ProcessGTExData.R_REQUIRED_FILES:
            if not os.path.exists( self.outputFolder + filePath):
                raise RainetException( "run_statistics : Input file is not present: " + filePath )
                
        # launch the analysis
        command = "cd " + os.path.dirname(ProcessGTExData.SWEAVE_R_SCRIPT) + \
                  "; Rscript %s %s %s %s %s %s %s %s" % ( 
                                                         ProcessGTExData.SWEAVE_R_SCRIPT, 
                                                         ProcessGTExData.WORKING_DIR,
                                                         self.outputFolder,
                                                         self.outputFolder+ProcessGTExData.ANNOTATION_OUTPUT_FILE, 
                                                         self.outputFolder+ProcessGTExData.EXPRESSION_OUTPUT_FILE,
                                                         self.outputFolder+ProcessGTExData.TX_EXPRESSION_OUTPUT_FOLDER,
                                                         self.outputFolder+ProcessGTExData.TX_EXPRESSION_AVG_OUTPUT_FILE,
                                                         self.outputFolder+ProcessGTExData.TX_EXPRESSION_AVG_OUTPUT_FILE_NO_OUTLIERS
                                                         )

        SubprocessUtil.run_command( command)


    # #
    # Method to remove files that are recreated each time.
    def clean_files(self):
        
        for filePath in ProcessGTExData.R_REQUIRED_FILES:
            files = glob.glob( self.outputFolder + filePath)
            for fi in files:
                os.remove(fi)


if __name__ == "__main__":
    
    try:
        
        #===============================================================================
        # Get input arguments, initialise class
        #===============================================================================
        parser = argparse.ArgumentParser(description='Script to produce basic statistics on the expression levels of transcripts, \
                                         normalise and modify the data (e.g. reduce sample dimension) for Rainet database insertion. ') 

        # positional args
        parser.add_argument('annotationFile', metavar='annotationFile', type=str,
                             help='GTEx file containing sample ID to tissue name mapping. E.g. GTEx_Data_V6_Annotations_SampleAttributesDS.txt')
        parser.add_argument('expressionFile', metavar='expressionFile', type=str,
                             help='GTEx file containing expression data per transcript per sample. E.g. GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt')
        parser.add_argument('tissueLevel', metavar='tissueLevel', type=str, choices=['SMTSD','SMTS'], 
                             help='Level of tissue / body-parts detail to be used. E.g. use of Brain, or of subtypes of brain.')
        parser.add_argument('minimumSamples', metavar='minimumSamples', type=int,
                             help='Minimum number of samples per tissue for tissue to be included')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Folder where to write output files.')
        # optional args
        parser.add_argument('--writeReport', metavar='writeReport', default = 1, type=int, help='Whether to write a pdf report using sampled data. (1 = Yes, 0 = No; Default = 1)')
        
        #display help when misusage
        if len(sys.argv) < 6: 
            parser.print_help()
    
        #gets the arguments
        args = parser.parse_args() 

        # add / to end of outputFolder
        if not args.outputFolder.endswith("/"):
            args.outputFolder+= "/"

        # Initialise class
        run = ProcessGTExData( args.annotationFile, args.expressionFile, args.tissueLevel, args.minimumSamples, args.outputFolder, args.writeReport )

        #===============================================================================
        # Run analysis / processing
        #===============================================================================
         
        # Create Logger instance by using the first log action.
        Logger.get_instance().info( "ProcessGTExData : Starting..." )
 
        # Start chrono
        Timer.get_instance().start_chrono()
 
        # Wipe output and log from previous runs
        run.clean_files()
 
        # Read annotation file
        Timer.get_instance().step( "reading annotation file..")    
        sampleTissue, tissueSample, problematicSamples = run.read_tissue_annotations()
        
        Timer.get_instance().step( "reading expression file..")    
        # Read expression file, using annotations and Process sample data into a single value
        # Note: nested function so that resource-heavy dictionary object is not duplicated
        run.average_sample_values(run.read_transcript_expression( sampleTissue, tissueSample, problematicSamples) )
         
        # Run R scripts / report
        run.run_statistics()


    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of ProcessGTExData. Aborting :\n" + rainet.to_string())

    # Stop the chrono      
    Timer.get_instance().stop_chrono( "ProcessGTExData : Finished" )

