package gulko;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.Scanner;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/*
 * Started 21-Feb-2017
 * @author Diogo Ribeiro
 * Script to process fitCons data from Gulko et al. 2015 Nat Genetics.
 * 
 * 1) Read gene/exon/CDS models in bed format, for which we want to know fitCons scores. Same features of the same transcript are merged (e.g. all CDS of transcript are merged together).
 * 2) Read fitCons scores from Gulko 2015 dataset, only retaining coordinates that are on the gene models we want.
 * 3) Produce output file with a line per transcript/feature, with number of nucleotides above a certain score (provided)
 */
public class ReadGulko {

	// global variables
	public static HashMap<String, String> coordTxDict;
	public static HashMap<String, ArrayList<String>> txCoordDict;
	public static HashMap<String, Double> coordScoreDict;

	// constructor
	public static void main(String[] args) {

		System.out.println("Starting..");

		// create Options object
		Options options = new Options();
		// add Gene models file argument
		options.addOption("geneModelsFile", true, "BED file containing gene models we want to get fitCons information for.");		
		options.addOption("gulkoFile", true, "Gulko 2015 file in bed format. Note: same genome assembly than geneModelsFile.");		
		options.addOption("outputFile", true, "Path of output file.");		
		
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = null;
		try { cmd = parser.parse( options, args);
		} catch (ParseException e) { e.printStackTrace(); }	
		
		// get options values
		String geneModelsFile = cmd.getOptionValue("geneModelsFile");
		String gulkoFile = cmd.getOptionValue("gulkoFile");
		String outputFile = cmd.getOptionValue("outputFile");

		// Input arguments
//		geneModelsFile = "/home/diogo/Documents/RAINET_data/lncRNA_info/gulko2013_fitCons/test/Rtest/all_exons_lincRNA_plus_CDS_head.bed";
//		gulkoFile = "/home/diogo/Documents/RAINET_data/lncRNA_info/gulko2013_fitCons/hg19_to_hg38/fc-i6-0_hg38.bed";

//		// testing
//		geneModelsFile = "/home/diogo/Documents/RAINET_data/lncRNA_info/gulko2013_fitCons/test/sampled_exons";
//		gulkoFile = "/home/diogo/Documents/RAINET_data/lncRNA_info/gulko2013_fitCons/test/fc-i6-0_shuf10000_sample.bed";		
		
		// Range of score thresholds to be used for analysis. Taken from Gulko et al. 2015 publication.
		double [] scoreRange = { 0.06, 0.07, 0.08, 0.09, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80};

		ReadGulko readGulko = new ReadGulko();

		System.out.println("Reading gene models..");
		readGulko.readGeneModelsFile( geneModelsFile);
		System.out.println("Reading Gulko2015 dataset..");
		readGulko.readGulkoFile( gulkoFile);
		System.out.println("Calculating scores per transcript..");
		readGulko.analyseTranscript( scoreRange, outputFile);

		System.out.println("FINISHED!");

	}

	
	/*
	 * Method to read file with gene models.
	 * 
	 * Produces dictionaries connecting transcripts to coordinates (1 nucleotide-based).
	 * If a transcript has several exons or other features, these are merged under the same transcript ID. 
	 * 
	 * @param geneModelsFile : bed file containing gene models we want to get fitCons information for.
	 * 
	 */
	public void readGeneModelsFile( String geneModelsFile){

		// variable declarations
		String line, tag, txEntry, txID, chrom;
		int start,end;
		Scanner data = null;
		String[] spl = null;

		// Initialisation
		coordTxDict = new HashMap<String,String>();  // key -> coordinate, value -> txID
		txCoordDict = new HashMap<String, ArrayList<String>>(); // key -> txID, value -> array of coordinates        

		try {
			data = new Scanner(new File( geneModelsFile));
		}
		catch ( IOException e) {
			e.printStackTrace();
			System.exit(0);
		}

		// Read file
		// example
		// chr1    12178   12227   exon:ENST00000450305.2:2        .       +       HAVANA  exon    .       ID=exon:ENST00000450305.2:2;Parent=ENST00000450305.2;gene_id=ENSG00000223972.5;transcript_id=ENST00000450305.2;gene_type=transcribed_unprocessed_pseudogene;gene_status=KNOWN;gene_name=DDX11L1;transcript_type=transcribed_unprocessed_pseudogene;transcript_status=KNOWN;transcript_name=DDX11L1-001;exon_number=2;exon_id=ENSE00001671638.2;level=2;transcript_support_level=NA;havana_gene=OTTHUMG00000000961.2;havana_transcript=OTTHUMT00000002844.2;ont=PGO:0000005,PGO:0000019;tag=basic

		while(data.hasNext()) {
			line = data.nextLine();
			spl = line.split("\t");
			chrom = spl[0];
			start = Integer.parseInt( spl[1]);
			end = Integer.parseInt( spl[2]);
			//            tag = chrom+":"+start+"-"+end;

//			txEntry = spl[3].split(":")[1];
//			txID = txEntry.split("\\.")[0];

			txID = spl[3].split("\\.")[0];
			
			// only initialise list of coordinates if transcript has not appeared
			// exons are joined in same transcript
			if (txCoordDict.get(txID) == null){
				txCoordDict.put(txID, new ArrayList<String>());
			}


			// TODO: check well start/end coordinates of transcripts
			// For each nucleotide in the sequence range
			for (int i = start; i < end; i++){
				String newTag = chrom + ":" + i;

				coordTxDict.put(newTag, txID); 

				txCoordDict.get(txID).add(newTag); 
			}

		}

		data.close();

		System.out.println( "Number of coordinates (in gene models) read:\t" + coordTxDict.size());
		System.out.println( "Number of transcripts read:\t" + txCoordDict.size());

		//		for (Entry<String, ArrayList<String>> entry : txCoordDict.entrySet()){
		//			System.out.println(entry.getKey() + "/" + entry.getValue().size());
		//		}

		//        for (Entry<String, String> entry : coordTxDict.entrySet()){
		//            System.out.println(entry.getKey() + "/" + entry.getValue());
		//        }

	}


	/*
	 * Method to read file with fitCons score per coordinate from Gulko 2015.
	 * 
	 * Reads scores for each coordinate, only storing coordinates that are within a gene model read in readGeneModelsFile method.
	 * 
	 * @param Gulko 2015 file in bed format.
	 */
	public void readGulkoFile( String gulkoFile){

		// variable declarations
		String line, chrom, tag;
		Scanner data = null;
		int start,end;
		double score;
		String[] spl = null;

		coordScoreDict = new HashMap<String,Double>();  // key -> coordinate, value -> score

		try {
			data = new Scanner(new File( gulkoFile));
		}
		catch ( IOException e) {
			e.printStackTrace();
			System.exit(0);
		}

		// Read file
		//		chr9    69098853        69098860        0.162113
		//		chr4    116064093       116064166       0.294411
		// These are range of coordinates with the same score

		int count = 0;

		while(data.hasNext()) {
			line = data.nextLine();
			spl = line.split("\t");
			chrom = spl[0];
			start = Integer.parseInt( spl[1]);
			end = Integer.parseInt( spl[2]);

			count++;
			if (count % 100000 == 0){
				System.out.printf( "readGulkoFile: Processed %s lines.. %n", count);
			}

			for (Integer i = start; i < end; i++){
				tag = chrom + ":" + i;

				// store only scores if coordinate exist in our transcripts
				if (coordTxDict.get( tag) != null){
					score = Double.parseDouble( spl[3]);
					coordScoreDict.put(tag, score); 
				}

			}
		}

		System.out.println( "Number of coordinates with score (in gene models):\t" + coordScoreDict.size());
		data.close();
		
	} // end readGulkoFile

	/*
	 * Method to produce statistics on fitCons score for each transcript.
	 * 
	 * For each fitCons score threshold, for each transcript, calculate how many nucleotides are above the score threshold.
	 * Note: if a coordinate belonging to the transcript does not contain any fitCons score information, that coordinate is not considered for any metric.
	 * 
	 * @param range of fitCons score thresholds
	 */
	public void analyseTranscript( double[] scoreRange, String outputFile){

		// create output file
		PrintWriter writer = null;
		try{
			writer = new PrintWriter(outputFile, "UTF-8"); 
			writer.println("score_cutoff\ttranscript\tncAbove\tncBelow\ttxSize\tpercAbove");
		} catch ( IOException e) {
			e.printStackTrace();
			System.exit(0);
		}

		// for each score threshold
		for ( double sco : scoreRange){

			// for each transcript
			for(Entry<String, ArrayList<String>> tx : txCoordDict.entrySet()) {
				String txID = tx.getKey();

				// count number of coordinates above and below/equal threshold
				int ncAbove = 0;
				int ncBelow = 0;
				// get transcript size
				int txSize = tx.getValue().size();

				// for each coordinate in transcript
				for (String coord : tx.getValue()){
					// if coordinate has any score
					if (coordScoreDict.get( coord) != null){
						// if score passes the current threshold
						if ( coordScoreDict.get( coord) > sco){ ncAbove++;}
						else{ ncBelow++;}
					} else{
						// some coordinates have no associated score. transcript size is reduced accordingly
						txSize--;
					}
				}


				if ( txSize != ( ncAbove + ncBelow)){
					System.out.println( "Wrong calculation of percentage above:\t" + txID);
					System.exit(1);
				}

				double percAbove = (double) ncAbove / (double) txSize;	

				//System.out.printf( "%s\t%s\t%s\t%s\t%s\t%.2f%n", sco, txID, ncAbove, ncBelow, txSize, percAbove);
				writer.printf("%s\t%s\t%s\t%s\t%s\t%.2f%n", sco, txID, ncAbove, ncBelow, txSize, percAbove);


			}

		}

		writer.close();

		// Validation
		//		0.06	ENST00000621762	302	56	358	0.84
		//
		//		grep ENST00000621762 exons_chr22 
		//		chr22	18501793	18501849	exon:ENST00000621762.1:2	.	-	HAVANA	exon	.	ID=exon:ENST00000621762.1:2;Parent=ENST00000621762.1;gene_id=ENSG00000278044.1;transcript_id=ENST00000621762.1;gene_type=lincRNA;gene_status=KNOWN;gene_name=XXbac-B33L19.12;transcript_type=lincRNA;transcript_status=KNOWN;transcript_name=XXbac-B33L19.12-001;exon_number=2;exon_id=ENSE00003730879.1;level=2;transcript_support_level=3;havana_gene=OTTHUMG00000188351.1;havana_transcript=OTTHUMT00000476955.1;tag=basic
		//		chr22	18511885	18512187	exon:ENST00000621762.1:1	.	-	HAVANA	exon	.	ID=exon:ENST00000621762.1:1;Parent=ENST00000621762.1;gene_id=ENSG00000278044.1;transcript_id=ENST00000621762.1;gene_type=lincRNA;gene_status=KNOWN;gene_name=XXbac-B33L19.12;transcript_type=lincRNA;transcript_status=KNOWN;transcript_name=XXbac-B33L19.12-001;exon_number=1;exon_id=ENSE00003725109.1;level=2;transcript_support_level=3;havana_gene=OTTHUMG00000188351.1;havana_transcript=OTTHUMT00000476955.1;tag=basic
		//
		//		chr22	18500417	18502257	0.053691
		//		chr22	18511600	18512000	0.095506
		//		chr22	18512000	18512080	0.119485
		//		chr22	18512080	18512098	0.077194
		//		chr22	18512098	18512188	0.116174
		//
		//		Summ:
		//		- 56 positions with 0.053691
		//		- 115 positions with 0.095506
		//		- 80 positions with 0.119485
		//		- 18 positions with 0.077194
		//		- 89 positions with 0.116174


	} // analyseTranscript end


} // class end

