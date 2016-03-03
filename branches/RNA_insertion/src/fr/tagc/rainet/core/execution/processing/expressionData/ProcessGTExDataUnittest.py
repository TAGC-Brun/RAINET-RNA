
import unittest

from ProcessGTExData import ProcessGTExData

# #
# Unittesting ProcessGTExData script
#
class ProcessGTExDataUnittest(unittest.TestCase):

    # #    
    # Run before each test
    def setUp(self):

        annotationFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/input_data/RNA/GTEx/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
        expressionFile = "/home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/testing_input_data/RNA/expression_test/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm_SHUFFLED.txt"
        tissueLevel = "SMTSD"
        minimumSamples = 30
        outputFolder = "unittest_folder/"

        self.run = ProcessGTExData(annotationFile, expressionFile, tissueLevel, minimumSamples, outputFolder, 0)

    # #
    # Test correct processing of GTEx_Data_V6_Annotations_SampleAttributesDS.txt file
    def test_annotation_parsing(self):
        
        sampleTissue, tissueSample, problematicSamples = self.run.read_tissue_annotations()
                
        # Some items are not real samples, or are missing information
        self.assertTrue( "NA12891-SM-5DWR9" in problematicSamples)

        # wc -l GTEx_Data_V6_Annotations_SampleAttributesDS.txt 
        # 11984 GTEx_Data_V6_Annotations_SampleAttributesDS.txt
        # Note: this includes header
        self.assertTrue( len(sampleTissue) == 11983 - len( problematicSamples))
        
        # grep GTEX-X88G-0426-SM-47JZ5 GTEx_Data_V6_Annotations_SampleAttributesDS.txt 
        # GTEX-X88G-0426-SM-47JZ5    0    B1    2 pieces  ~5mm d.  Adherent fibroadipose tissue, discontinuous, delineated, up to ~2mm    6.8    Blood Vessel    Artery - Tibial    0007610    1139    Presumed Death    BP-32413    RNA isolation_PAXgene Tissue miRNA    10/22/2012    LCSET-2974    USE ME        0.93951505    111923    0.97698283    854    0.94767725    0.8776571    0.28182453    24494    0.7805154    0.7434379    76    54.057327    0.0026057002    211    69254767    68577653    0.022762924    378318    72363938    7905460    0.63382757    142291    33274718    57811    0.032747243    0    0.09932569    0.73967665    0.8317357    53525915    0.0022174604    14127206    2830916    14214546    202    14357831    12130061    0.0024132524    14476618    50.205982    0.005227991    0.95583934    837    0.21948458    50.154083        
        self.assertTrue( sampleTissue["GTEX-X88G-0426-SM-47JZ5"] == "Artery - Tibial")

        # cut -f7 GTEx_Data_V6_Annotations_SampleAttributesDS.txt | sort -u | wc -l
        # 56
        self.assertTrue( len(tissueSample) == 56 -2) # -2 refers to the header and an exceptional blank line
        
        # cut -f7 GTEx_Data_V6_Annotations_SampleAttributesDS.txt | grep "Thyroid" | wc -l
        # 437
        self.assertTrue( len(tissueSample["Thyroid"]) == 437)


    # #
    # Test correct processing of GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt file
    def test_expression_parsing(self):

        sampleTissue, tissueSample, problematicSamples = self.run.read_tissue_annotations()
        
        txExpressionTissue = self.run.read_transcript_expression( sampleTissue, tissueSample, problematicSamples) 

        # wc -l GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm_SHUFFLED.txt 
        # 51 GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm_SHUFFLED.txt
        self.assertTrue( len(txExpressionTissue) == 51)

        # Note: there is 49 tissue that are kept with the used parameters. 
        # one of the 54 previous tissues simply do not show up in the expression file ?   
        # 4 of the remaining 53 tissue are excluded with minimumSample cutoff == 30 
        self.assertTrue( len(txExpressionTissue["ENST00000461495"]) == 49 ) 

        allValues = set()
        for tiss in txExpressionTissue["ENST00000461495"]:
            for val in txExpressionTissue["ENST00000461495"][tiss]:
                allValues.add( val)

        # with grep, looked at a non-zero value                
        self.assertTrue( 0.886113 in allValues ) 


    # #
    # Test IQR filtering function
    def test_iqr_filtering(self):
        
        fakeArray = [0,0,0,10,10,10,100,100]
    
        valuesInRange = self.run.remove_outliers( fakeArray)

        # Checked if same result as in R
        # > fakeArray <- c( 0,0,0,10,10,10,100,100)
        # > 
        # > IQR(fakeArray)
        # [1] 32.5
        # > 
        # > quantile(fakeArray)
        #    0%   25%   50%   75%  100% 
        #   0.0   0.0  10.0  32.5 100.0 
        # > 
        # > q1 = 0.0
        # > q3 = 32.5
        # > 
        # > q1 - 1.5*IQR(fakeArray)
        # [1] -48.75
        # > q3 + 1.5*IQR(fakeArray)
        # [1] 81.25        

        self.assertTrue( 100 not in valuesInRange)


    # #
    # Test correct production of output files to be used for insertion
    # Note: I will not test files that are used for report, which are randomised each run
    def test_average_sample_values(self):
        
        sampleTissue, tissueSample, problematicSamples = self.run.read_tissue_annotations()        
        txExpressionTissue = self.run.read_transcript_expression( sampleTissue, tissueSample, problematicSamples) 
        
        self.run.average_sample_values(txExpressionTissue)

        # for simplicity, I confirm that the values before outlier filtering are correct. Since outlier filter has been confirm independently
        # Using R
        # > realArray <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.65244899999999995, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.64329599999999998, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0922700000000001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        # > mean(realArray)
        # [1] 0.01388381
        # > sd <- sqrt(sum((realArray - mean(realArray))^2) / (length(realArray) ))
        # > sd
        # [1] 0.1078171
        # > median(realArray)
        # [1] 0
        # > sd/mean(realArray)
        # [1] 7.765672
        # > max(realArray)
        # [1] 1.09227

        # Note: in python numpy, standard deviation is calculated by sqrt(sum((x - mean(x))^2) / N (i.e. not N-1) 
       
        validatedValues = ["%.3f" % 0.01388381, "%.3f" %  0.1078171, "%.3f" % 0, 7.765672, "%.3f" % 1.09227]
        
        with open(self.run.outputFolder+"/transcript_expression_metrics.tsv") as f:
            for line in f:
                if "ENST00000461495" in line and "Testis" in line:
                    #tx, tiss, mean, sd, median, coef, max 
                    spl = line.split("\t")[2:]

                    for i in xrange(len(spl)):               
                        
                        self.assertTrue("%.3f" % float(spl[i])  == "%.3f" % float(validatedValues[i]) )


    # #    
    # Run after each test    
    def tearDown(self):
        pass
    

