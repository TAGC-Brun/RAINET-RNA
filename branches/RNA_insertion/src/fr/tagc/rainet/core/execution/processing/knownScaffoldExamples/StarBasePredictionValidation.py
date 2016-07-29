import sys
import os
import argparse
import glob
import numpy as np
import random 
import pandas as pd
from scipy import stats

from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil

from sqlalchemy import or_,and_
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinInteraction import ProteinInteraction
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.RNACrossReference import RNACrossReference
from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.util.data.DataManager import DataManager
from NPInterPredictionValidation import NPInterPredictionValidation

#===============================================================================
# Started 20-Apr-2016 
# Diogo Ribeiro
# Script to see if catRAPID predictions distinguish starBase CLIP interactions. Based on NPInterPredictionValidation.py
#===============================================================================

#===============================================================================
# General plan:
# 1) Read starBase file, apply CLIP read number filters
# 2) Use uniprot mapping of starbase protein names. Use RAINET database (currently Ensembl v82) to map transcript names to Ensembl IDs
# 3) Read catRAPID omics file, store pairs and their interaction scores
# 4) Use R for plotting distributions
#===============================================================================
#===============================================================================
# Processing notes:
# 
# - StarBase provides data on RBP gene symbol and LncRNA gene symbol, not on transcript level (even though they map against transcript).
#   to process their data, we assume that if a CLIP read maps to a transcript of a gene, it will map as well for all transcripts of that gene (which may not be true)
# - In output files, if same RNA-protein pair is repeated (due to input or mapping), the higher clipNumber and BioComplex values are kept.
#===============================================================================


class StarBasePredictionValidation( object ):

    DISTRIBUTION_SCRIPT = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/knownScaffoldExamples/StarBase_stats.R"

    def __init__(self, catrapidFile, starbaseFile, starbaseConversionFile, rainetDB, outputFolder, minimumBioComplex, minimumClipReadNumber):

        self.catrapidFile = catrapidFile
        self.starbaseFile = starbaseFile
        self.starbaseConversionFile = starbaseConversionFile
        self.rainetDB = rainetDB
        self.outputFolder = outputFolder
        self.minimumBioComplex =  minimumBioComplex
        self.minimumClipReadNumber = minimumClipReadNumber
    
        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()

        # make output folder
        if not os.path.exists( self.outputFolder):
            os.mkdir( self.outputFolder)


    # #
    # Use RAINET DB to retrieve Protein cross references
    def protein_cross_references(self):
        
        # # Get external references
        # produce dictionary where key is xref ID and value the uniprotAC
        # Note: an external ID can point to several uniprot IDs        
        query = self.sql_session.query( ProteinCrossReference.protein_id, ProteinCrossReference.crossReferenceID ).filter( ProteinCrossReference.sourceDB == "Ensembl_PRO" ).all()

        protCrossReference = {} # key -> external ID, val -> set of uniprotACs        
        for uniprotID, externalID in query:     
            if externalID not in protCrossReference:
                protCrossReference[ externalID] = set()
            protCrossReference[ externalID].add( str( uniprotID))

        return protCrossReference


    # #
    # Use RAINET DB to retrieve RNA cross references
    def rna_cross_references(self):
        
        query = self.sql_session.query( RNA.transcriptID, RNA.externalGeneName ).all()
        
        rnaCrossReference = {} # key -> gene name, val -> list of ensembl IDs
        # Note: an gene name points to several ensembl IDs        
        
        for ensemblID, geneName in query:
            if geneName not in rnaCrossReference:
                rnaCrossReference[ geneName] = []
                
            rnaCrossReference[ geneName].append( str( ensemblID) )

        return rnaCrossReference


    # #
    # Read uniprot file that maps starbase gene names into uniprotAC
    def read_starbase_protein_conversion_file(self):
        
        # Header: yourlist:M2016042065ZNJ43KQK    Entry   Entry name      Status  Protein names   Gene names      Organism        Length


        inFile = open( self.starbaseConversionFile, "r")

        conversionDict = {} # key -> starbase gene name, value -> list of uniprotACs
        for line in inFile:
            line = line.strip()
            spl = line.split( "\t")

            originalIDs = spl[0].split( ",")
            uniprotID = spl[1]

            # Note: sometimes a single line includes two items e.g.:
            # LIN28,LIN28A    Q9H9Z2  LN28A_HUMAN     reviewed        Protein lin-28 homolog A (Lin-28A) (Zinc finger CCHC domain-containing protein 1)       LIN28A CSDD1 LIN28 ZCCHC1       Homo sapiens (Human)    209            
            for originalID in originalIDs:
                            
                if originalID not in conversionDict:
                    conversionDict[ originalID] = []
                
                conversionDict[ originalID].append( uniprotID)
        
        return conversionDict

        
    # #
    # Read starbase file and retrieve list of proteins interacting with wanted RNA
    def read_starbase_file(self):
        
        # Note: StarBase uses Gene Symbols for Proteins, and Gene Symbols for transcripts
        # Header : name    geneName        targetSites     bioComplex      clipReadNum     CancerNum

        #===================================================================
        # Read Starbase file using header line 
        #===================================================================
        table = pd.read_table( self.starbaseFile, header = 0, sep = "\t", skip_blank_lines = True)
 
        print "read_starbase_file: Number interactions before any filter:",len(table)
  
        filteredTable = table.copy()
          
        #===================================================================
        # Field filtering on StarBase data
        #===================================================================
          
        filteredTable = filteredTable.loc[filteredTable["bioComplex"] >= self.minimumBioComplex]
           
        filteredTable = filteredTable.loc[filteredTable["clipReadNum"] >= self.minimumClipReadNumber] 
  
        print "read_starbase_file: Number interactions after bioComplex and clipReadNumber filter:",len(filteredTable)
 
        #===================================================================
        # Retrieve interactions
        #===================================================================        
        
        missedProteins = set()
        missedRNAs = set()
 
        interactingPairs = {} # key -> pair of transcriptID and proteinID, val -> number of CLIP reads mapped
        interactingPairsBiocomplex = {} # key -> pair of transcriptID and proteinID, val -> number of biocomplexes
        setOfRNAs = set()
        setOfProts = set()
 
        for index, row in filteredTable.iterrows():
            starbaseProtID = str( row[ "name"]).upper()
            starbaseRNAID = str( row[ "geneName"]).upper()
            
            clipReads = int( row["clipReadNum"])
            bioComplex = int( row["bioComplex"])

            # dealing with naming exceptions            
            if starbaseProtID == "FMRP":
                starbaseProtID = "FMR1"
            if starbaseProtID == "EIF4AIII":
                starbaseProtID = "EIF4A3"

            # checking if ID can be traced
            boo = 1            
            if starbaseProtID not in self.starbaseProtXrefDict:
                missedProteins.add( starbaseProtID)
                boo = 0
            if starbaseRNAID not in self.starbaseRNAXrefDict:
                missedRNAs.add( starbaseRNAID)
                boo = 0
            if boo == 0:
                continue

            # normally there is several IDs both for Gene name and uniprotAC, and for Gene name and transcripts  
            # we expand the interaction for all possible ID combinations     
            uniprotIDs = self.starbaseProtXrefDict[ starbaseProtID]
            ensemblIDs = self.starbaseRNAXrefDict[ starbaseRNAID]

            for uniprotID in uniprotIDs:
                for ensemblID in ensemblIDs:
                    tag = ensemblID + "|" + uniprotID
                    
                    if tag not in interactingPairs:
                        interactingPairs[ tag] = 0
                        interactingPairsBiocomplex[ tag] = 0
                    # keep the maximum value of clip reads mapped
                    if clipReads > interactingPairs[ tag]:
                        interactingPairs[ tag] = clipReads
                    if bioComplex > interactingPairsBiocomplex[ tag]:
                        interactingPairsBiocomplex[ tag] = bioComplex
        
                    setOfProts.add( uniprotID)
                    setOfRNAs.add( ensemblID)
        
        print "Missing proteins:",len(missedProteins), missedProteins
        print "Missing RNAs",len(missedRNAs), missedRNAs

        assert len( interactingPairs) == len( interactingPairsBiocomplex)

        print "read_starbase_file: Total number of interacting pairs:",len(interactingPairs)
        print "read_starbase_file: Total number of interacting RNAs:",len(setOfRNAs)
        print "read_starbase_file: Total number of interacting proteins:",len(setOfProts)
      
      
        self.starbasePairs = interactingPairs
        return interactingPairs, interactingPairsBiocomplex


    # #
    # Read catRAPID file, match peptide IDs to protein IDs, retrieve scores.
    def read_catrapid_file(self):


        interactingPairs = {} # key -> pair of transcriptID and proteinID, val -> maximum interaction score
        
        #e.g. 1       1       ENSP00000269701_ENST00000456726 -266.23 0.986

        peptideIDNotFound = set()
        proteinSet = set()

        countLines = 0
        
        with open( self.catrapidFile, "r") as f:
            for line in f:
                spl = line.split( "\t")

                countLines+= 1 

                if countLines % 10000000 == 0:
                    print "Processed %s interactions" % countLines
                               
                splIDs = spl[2].split("_")
                peptideID = splIDs[0]
                transcriptID = splIDs[1]
                intScore = float( spl[3] )

                if peptideID in self.xrefDict:
                    proteinID = self.xrefDict[ peptideID]
                    if len( proteinID) == 1:
                        #proteinID = next( iter( proteinID))
                        proteinID, = proteinID # unpacking set
                    else:
                        raise RainetException( "ENSP should point to a single UniProtID: " + line)                        
                else:
                    #print "read_catrapid_file: PeptideID not found in RAINET database: ", peptideID 
                    peptideIDNotFound.add( peptideID)
                    continue

                pair = transcriptID + "|" + proteinID

                proteinSet.add(proteinID)

                # add pair to interacting pairs and keep the maximum interaction score
                if pair not in interactingPairs:
                    interactingPairs[ pair] = float("-inf")
                if intScore > interactingPairs[ pair]:
                    interactingPairs[ pair] = intScore

        print "read_catrapid_file: Number of peptideIDs not found in RAINET DB: ", len( peptideIDNotFound) # for old catRAPID dataset, 243 is expected
        print "read_catrapid_file: Number of proteins: ", len( proteinSet)
        print "read_catrapid_file: Number of protein-RNA pairs in catRAPID: ", len( interactingPairs)

        return interactingPairs


    # #
    # Read catRAPID file, match peptide IDs to protein IDs, retrieve scores.
    # Version to read a slightly different format, used for the mRNA dataset
    def read_catrapid_file_mRNA(self):

        # E.g.: ENSP00000001008    ENST00000000233    -62.65    0.00    0.09    0.17    0.11    0.14

        interactingPairs = {} # key -> pair of transcriptID and proteinID, val -> maximum interaction score
        
        peptideIDNotFound = set()
        proteinSet = set()
        
        with open( self.catrapidFile, "r") as f:
            for line in f:
                spl = line.split( "\t")
                
                peptideID = spl[0]
                transcriptID = spl[1]
                intScore = float( spl[2] )

                if peptideID in self.xrefDict:
                    proteinID = self.xrefDict[ peptideID]
                    if len( proteinID) == 1:
                        #proteinID = next( iter( proteinID))
                        proteinID, = proteinID # unpacking set
                    else:
                        raise RainetException( "ENSP should point to a single UniProtID: " + line)                        
                else:
                    #print "read_catrapid_file: PeptideID not found in RAINET database: ", peptideID 
                    peptideIDNotFound.add( peptideID)
                    continue

                pair = transcriptID + "|" + proteinID

                proteinSet.add(proteinID)

                # add pair to interacting pairs and keep the maximum interaction score
                if pair not in interactingPairs:
                    interactingPairs[ pair] = float("-inf")
                if intScore > interactingPairs[ pair]:
                    interactingPairs[ pair] = intScore

        print "read_catrapid_file: Number of peptideIDs not found in RAINET DB: ", len( peptideIDNotFound)
        print "read_catrapid_file: Number of proteins: ", len( proteinSet)
        print "read_catrapid_file: Number of protein-RNA pairs in catRAPID: ", len( interactingPairs)

        return interactingPairs


    # #
    # Read catRAPID file.
    # Updated for new catRAPID format. No need for cross references.
    def read_catrapid_file_new(self):

        # E.g.: sp|Q6P6C2|ALKB5_HUMAN ENST00000559683   47.85   0.93    0.23

        interactingPairs = {} # key -> pair of transcriptID and proteinID, val -> score
        
        proteinSet = set()

        countLines = 0
        
        with open( self.catrapidFile, "r") as f:
            for line in f:
                spl = line.split(" ")

                countLines+= 1 

                if countLines % 10000000 == 0:
                    print "Processed %s interactions" % countLines
                               
                
                proteinID = spl[0].split( "|")[1]
                spl2 = spl[1].split( "\t")
                transcriptID = spl2[0]
                intScore = float( spl2[1])
                
                pair = transcriptID + "|" + proteinID

                proteinSet.add(proteinID)

                # add pair to interacting pairs and keep the maximum interaction score
                if pair not in interactingPairs:
                    interactingPairs[ pair] = intScore
                else:
                    raise RainetException( "Repeated protein-RNA pair: " + line)
   

        print "read_catrapid_file_new: Number of proteins: ", len( proteinSet)
        print "read_catrapid_file_new: Number of protein-RNA pairs in catRAPID: ", len( interactingPairs)

        return interactingPairs


    # #
    # Read expression file instead of catrapid file.
    def read_expression_file(self):

        # E.g.: Q7RTM1  ENST00000437598 0.013

        interactingPairs = {} # key -> pair of transcriptID and proteinID, val -> score
        
        proteinSet = set()

        countLines = 0

        outFile = open( run.outputFolder + "/scores.tsv", "w")
        outFile.write( "pairID\tcatrapid_score\tin_validated_set\n")

        text = ""

        with open( self.catrapidFile, "r") as f:
            for line in f:
                spl = line.split("\t")

                countLines+= 1 

                if countLines % 1000000 == 0:
                    print "Processed %s interactions" % countLines
                    
                    outFile.write(text)
                    
                    interactingPairs = {}
                    text = ""
                
                proteinID = spl[0]
                transcriptID = spl[1]
                intScore = float( spl[2])
                
                pair = transcriptID + "|" + proteinID

                proteinSet.add(proteinID)

                # add pair to interacting pairs and keep the maximum interaction score
                if pair not in interactingPairs:
                    interactingPairs[ pair] = intScore
                else:
                    raise RainetException( "Repeated protein-RNA pair: " + line)

                if pair in self.starbasePairs:
                    inValidated = 1
                else:
                    inValidated = 0
   
                text += "%s\t%s\t%s\n" % ( pair, intScore, inValidated) 

        outFile.write(text)
        outFile.close()
   
        print "read_expression_file: Number of proteins: ", len( proteinSet)
        print "read_expression_file: Number of protein-RNA pairs in file: ", len( interactingPairs)

        return interactingPairs


if __name__ == "__main__":
    
    try:
        # Create Logger instance by using the first log action.
        Logger.get_instance().info( "StarBasePredictionValidation : Starting..." )

        #===============================================================================
        # Get input arguments, initialise class
        #===============================================================================
        parser = argparse.ArgumentParser(description='# Script to see if catRAPID predictions distinguish Starbase interactions ') 

        # positional args
        parser.add_argument('catRAPIDFile', metavar='catRAPIDFile', type=str,
                             help='File path of CatRAPID omics/fragments results from the webserver.')
        parser.add_argument('starBaseFile', metavar='starBaseFile', type=str,
                             help='File path of StarBase file.')
        parser.add_argument('starBaseProteinConversionFile', metavar='starBaseProteinConversionFile', type=str,
                             help='File path of StarBase file.')
        parser.add_argument('rainetDB', metavar='rainetDB', type=str, help='Path to RAINET database to be used.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Folder where to write output files.')

        parser.add_argument('--minimumBioComplex', metavar='minimumBioComplex', default = 0, type=int, help='Minimum value of BioComplex to keep starBase interaction.')
        parser.add_argument('--minimumClipReadNumber', metavar='minimumClipReadNumber', default = 0, type=int, help='Minimum value of clipReadNum to keep starBase interaction.')
        parser.add_argument('--catRAPIDmRNA', metavar='catRAPIDmRNA', default = 0, type=int, help='Whether provided catRAPID file is mRNA or lncRNA file. Parser used is different.')
        parser.add_argument('--newFormat', metavar='newFormat', default = 0, type=int, help='Whether provided catRAPID file is using new format (e.g. with uniprotac instead of ENSP) or not. If value == 2, read expression file instead.')

        
        #gets the arguments
        args = parser.parse_args( ) 

        # Initialise class
        run = StarBasePredictionValidation( args.catRAPIDFile, args.starBaseFile, args.starBaseProteinConversionFile, args.rainetDB, 
                                            args.outputFolder, args.minimumBioComplex, args.minimumClipReadNumber)

        #===============================================================================
        # Run analysis / processing
        #===============================================================================
         
        # Start chrono
        Timer.get_instance().start_chrono()
 
        Timer.get_instance().step( "reading RAINET DB and conversion files..")    

        # Build starbase RNA cross references
        run.starbaseRNAXrefDict = run.rna_cross_references()
 
        # Build starbase Protein cross references
        run.starbaseProtXrefDict = run.read_starbase_protein_conversion_file()
        
        # Build Protein cross references
        run.xrefDict = run.protein_cross_references()
 
        Timer.get_instance().step( "reading starBase file..")    

        starbasePairs, starbasePairsBiocomplex = run.read_starbase_file()

        # #
        # Write file with all Starbase pairs, regardless of being in catRAPID predictions or not
        outFile = open( run.outputFolder + "/starbase_interactions.tsv", "w")
        outFile.write( "pairID\tclipReads\tbioComplex\n")
        
        for pair in starbasePairs:
            outFile.write( "%s\t%s\t%s\n" % (pair, starbasePairs[ pair], starbasePairsBiocomplex[ pair]) )            
        outFile.close()

        Timer.get_instance().step( "reading catRAPID file..")    

        if args.newFormat == 1:
            catrapidPairs = run.read_catrapid_file_new()
        elif args.newFormat == 2:
            catrapidPairs = run.read_expression_file()
        elif args.catRAPIDmRNA:
            catrapidPairs = run.read_catrapid_file_mRNA()
        else: 
            catrapidPairs = run.read_catrapid_file()

        if args.newFormat != 2:

            # #
            # Quick stats on the data
            countYes = 0
            countNo = 0
            countYesSum = 0
            countNoSum = 0
            for pair in catrapidPairs:
                if pair in starbasePairs:
                    countYes+=1
                    countYesSum+= catrapidPairs[ pair]
                else:
                    countNo+=1
                    countNoSum+= catrapidPairs[ pair]
      
            print "True: %s\tFalse: %s" % ( countYes, countNo)
            print "True sum: %s\tFalse sum: %s" % ( countYesSum, countNoSum)
            print "True mean: %.2f\tFalse mean: %.2f" % ( countYesSum / float( countYes), countNoSum / float( countNo))
    
            Timer.get_instance().step( "Writing output file..")    
     
            # #
            # Write file for R processing
            outFile = open( run.outputFolder + "/scores.tsv", "w")
            outFile.write( "pairID\tcatrapid_score\tin_validated_set\tclipReads\n")
             
            # array of 1s and 0s, whether pair in Starbase or not
            inValidated = [ 1 if pair in starbasePairs else 0 for pair in catrapidPairs]
     
            for i, pair in enumerate( catrapidPairs):
                if pair in starbasePairs:
                    outFile.write( "%s\t%s\t%s\t%s\n" % ( pair, catrapidPairs[ pair], inValidated[ i], starbasePairs[ pair] ) )
                else:
                    outFile.write( "%s\t%s\t%s\t%s\n" % ( pair, catrapidPairs[ pair], inValidated[ i], 0 ) )
    
            outFile.close()
    
    
    #         # # Run R command to create figure
    #         command = "Rscript %s %s %s" % ( StarBasePredictionValidation.DISTRIBUTION_SCRIPT, outFile.name, run.outputFolder)
    #         result = SubprocessUtil.run_command( command) #, return_stdout = 1, verbose = 1)


    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of StarBasePredictionValidation. Aborting :\n" + rainet.to_string())

    # Stop the chrono      
    Timer.get_instance().stop_chrono( "StarBasePredictionValidation : Finished" )


