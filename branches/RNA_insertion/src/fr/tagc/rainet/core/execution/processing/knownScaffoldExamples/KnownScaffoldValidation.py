
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


#===============================================================================
# Started 21-Mar-2016 
# Diogo Ribeiro
# Script to see if known scaffolding RNAs have enriched binding to their known protein targets when using catRAPID data
#
#===============================================================================

#===============================================================================
# General plan:
# 
# To be run for each RNA of interest:
# - Read supposed interacting proteins from NPInter database file or manual list
# - Match those proteins to uniprotAC using RAINET database (if needed)
# - Read catRAPID omics/fragments file, apply interaction cutoff 
# - Perform enrichment analysis 
#===============================================================================

#===============================================================================
# Processing notes:
# 
# - If reading catRAPID fragments, we retain RNA-protein interaction if at least one of the fragments passes the score cutoffs
# - Assuming that NPInter file has RNA as moleculeA and partner as moleculeB (the case for all examples I picked so far)
# - Running various parameters at: /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/results/knownScaffoldValidation
#===============================================================================


class KnownScaffoldValidation( object ):

    SPECIES = "Homo sapiens"
    MOLECULEBTYPE = "protein"
    HYPERGEOMETRIC_TEST_SCRIPT = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/knownScaffoldExamples/hypergeometric_test.R"
    DISTRIBUTION_SCRIPT = "/home/diogo/workspace/tagc-rainet-RNA/src/fr/tagc/rainet/core/execution/processing/knownScaffoldExamples/plot_distribution.R"
    CATRAPID_ZS_MEAN = 23.25
    CATRAPID_ZS_STD = 37.90

    def __init__(self, catRAPIDFile, validatedFile, wantedRNAFile, rainetDB, outputFolder, 
                 discriminativePowerCutoff, zscoreCutoff, topProportion, npinter, columnForPlot, searchSpaceFile, zscoreToPropensity):

        self.catRAPIDFile = catRAPIDFile
        self.validatedFile = validatedFile
        self.wantedRNAFile = wantedRNAFile
        self.rainetDB = rainetDB
        self.outputFolder = outputFolder
        self.discriminativePowerCutoff = discriminativePowerCutoff
        self.zscoreCutoff = zscoreCutoff
        self.topProportion = topProportion
        self.npinter = npinter
        self.columnForPlot = columnForPlot
        self.searchSpaceFile = searchSpaceFile
        self.zscoreToPropensity = zscoreToPropensity

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()

        # Run query to retrieve all proteins (uniProtACs) in RAINET database
        self.proteinsInRainet, self.xrefDict = self._proteins_in_rainet()

    # #
    # Access RAINET database and retrieve all uniprotACs included in the database
    def _proteins_in_rainet(self):
        
        query = self.sql_session.query( Protein.uniprotAC ).all()
        
        proteinsInRainet = { str(prot[0]) for prot in query}

        # produce dictionary where key is xref ID and value the uniprotAC
        query = self.sql_session.query( ProteinCrossReference.protein_id, ProteinCrossReference.crossReferenceID ).all()
        xrefDict = { str(prot[1]) : str(prot[0]) for prot in query } # dict comprehension
       
        return proteinsInRainet, xrefDict

    # #
    # Read catRAPID file, apply wanted interaction score filters and retrieve set of interacting and non-interacting proteins
    def read_catRAPID_file(self):

        # Provided catRAPID header         # Note: catRAPID file header does not match content and is not tab-delimited
        ## Protein ID   Transcript ID   Z-score   Discriminative Power   Interaction Strength   Motif Reference   Domain Presence   Domain Score   -   Star Rating Score

        #===================================================================
        # read catRAPID file skipping header line 
        #===================================================================
        table = pd.read_table( self.catRAPIDFile, header = None, sep = "\t", skip_blank_lines = True, skiprows = 1)
        print "read_catRAPID_file: Starting number of interactions:",len(table)               
#         print table.head(5)
        
        #===================================================================
        # retrieve list of scores per target before any filtering, use maximum score between fragments
        #===================================================================
        
        scoreDict = {}
        for index, row in table.iterrows():
            # split different protein IDs and transcript
            spl = row[0].split("|")
            uniprotAC = spl[1]
            
            score = row[self.columnForPlot] 
        
            if uniprotAC not in scoreDict:
                scoreDict[uniprotAC] = float("-inf")
            # since there are a score for each RNA fragment, store only the maximum z-score among fragments
            if score > scoreDict[uniprotAC]:
                scoreDict[uniprotAC] = score

#         # new code for calculating average instead of max
#             if uniprotAC not in scoreDict:
#                 scoreDict[uniprotAC] = []
#             scoreDict[uniprotAC].append( score)
# 
#         for prot in scoreDict:
#             # calculate average
#             scoreDict[ prot] = sum( scoreDict[ prot]) / float( len( scoreDict[ prot]))

                
        #===================================================================
        # Filter for discriminative power
        #===================================================================
        filteredTable = table.loc[table.loc[:,2] > self.discriminativePowerCutoff,:]
        print "read_catRAPID_file: Number interactions after discriminative power filter:",len(filteredTable)

#         #===================================================================
#         # Filter for interaction strength
#         #===================================================================
#         filteredTable = filteredTable.loc[filteredTable.loc[:,3] > self.zscoreCutoff,:]
#         print "read_catRAPID_file: Number interactions after interaction strength filter:",len(filteredTable)
# 
#         # store names of catRAPID proteins not present in RAINET
#         missingProteins = set()
# 
#         headValue = int( "%i" % (len( filteredTable) * self.topProportion ) )
#         
#         # if the topProportion option was used, sort table and get only the top scores
#         if len(filteredTable) != headValue:
#         
#             # Sort first by interaction strength and then by discriminative power, higher to lower
#             sortedTable = filteredTable.sort_values([3, 2], ascending=[False, False])
#     
#             filteredTable = sortedTable.head(headValue)

        #===================================================================
        # Filter for Z-score 
        #===================================================================
        filteredTable = filteredTable.loc[filteredTable.loc[:,1] > self.zscoreCutoff,:]
        print "read_catRAPID_file: Number interactions after zscore filter:",len(filteredTable)

        # store names of catRAPID proteins not present in RAINET
        missingProteins = set()

        headValue = int( "%i" % (len( filteredTable) * self.topProportion ) )
        
        # if the topProportion option was used, sort table and get only the top scores
        if len(filteredTable) != headValue:
        
            # Sort first by zscore and then by discriminative power, higher to lower
            sortedTable = filteredTable.sort_values([1, 2], ascending=[False, False])
    
            filteredTable = sortedTable.head(headValue)

        #===================================================================
        # Get set of interacting proteins  
        #===================================================================
        # Note: the transcript may have been fragmented, so there will be several hits for same protein - transcript pair.

        interactingProts = {}

        # get item storing transcript-protein pair
        interactions = filteredTable.loc[:,0]
        
        for inter in interactions:
            # split different protein IDs and transcript
            spl = inter.split("|")
            uniprotAC = spl[1]
            
            if uniprotAC not in interactingProts:
                interactingProts[uniprotAC] = 0
            interactingProts[uniprotAC] += 1
            
            # If one wants to filter out proteins not in RAINET
#             if uniprotAC in self.proteinsInRainet:
#                 if uniprotAC not in interactingProts:
#                     interactingProts[uniprotAC] = 0
#                 interactingProts[uniprotAC] += 1
#             else:
#                 # this can be as proteinAC was deprecated etc
#                 # print "read_catRAPID_file: ProteinID not found in RAINET database", uniprotAC
#                 missingProteins.add( uniprotAC)

        #===================================================================
        # Get set of non-interacting proteins  
        #===================================================================
        # Note: loop over all data in catRAPID, check if in RAINET database, and keep only if not previously tagged as interacting

        setOfNonInteractingProts = {}

        # get item storing transcript-protein pair in original table
        interactions = table.loc[:,0]
        
        for inter in interactions:
            # split different protein IDs and transcript
            spl = inter.split("|")
            uniprotAC = spl[1]

            if uniprotAC not in interactingProts:
                if uniprotAC not in setOfNonInteractingProts:
                    setOfNonInteractingProts[uniprotAC] = 0                    
                setOfNonInteractingProts[uniprotAC] += 1

#             if uniprotAC in self.proteinsInRainet:
#                 # check if not tagged as interacting
#                 if uniprotAC not in interactingProts:
#                     if uniprotAC not in setOfNonInteractingProts:
#                         setOfNonInteractingProts[uniprotAC] = 0                    
#                     setOfNonInteractingProts[uniprotAC] += 1
#             else:
#                 # this can be as proteinAC was deprecated etc
#                 # print "read_catRAPID_file: ProteinID not found in RAINET database", uniprotAC
#                 missingProteins.add( uniprotAC)

        #print "read_catRAPID_file: Total number of proteins in catRAPID file not present in RAINET DB:",len(missingProteins)
        print "read_catRAPID_file: Total number of proteins in catRAPID file:",len(interactingProts)+len(setOfNonInteractingProts)
        print "read_catRAPID_file: Total number of interacting proteins:",len(interactingProts)
        print "read_catRAPID_file: Total number of non-interacting proteins:",len(setOfNonInteractingProts)
        
        return interactingProts, setOfNonInteractingProts, scoreDict


    # #
    # Read NPInter file and retrieve list of proteins interacting with wanted RNA
    def read_NPInter_file(self):
        
        # Note: NPInter uses NONCODE database for transcript IDs, but in fact uses their Gene IDs, not the transcript IDs
           
        #===================================================================
        # read wantedRNAFile,  file containing IDs to look for on NPInter file
        #===================================================================
        wantedRNATable = pd.read_table( self.wantedRNAFile, header = None, sep = "\t", skip_blank_lines = True)

        #===================================================================
        # read NPInter file using header line 
        #===================================================================
        table = pd.read_table( self.validatedFile, header = 0, sep = "\t", skip_blank_lines = True)

        print "read_NPInter_file: Number interactions before any filter:",len(table)
 
        #===================================================================
        # get lines of NPInter file containing wanted RNA in wanted column
        #===================================================================        
        wantedEntries = pd.DataFrame()
        
        for row in wantedRNATable.iterrows():
            key = row[1][0]
            val = row[1][1]
            entries = table.loc[table[key] == val]

            if len(wantedEntries) == 0:
                wantedEntries = entries.copy()
            else:
                # note that pandas append returns object and does not change object in place
                wantedEntries = wantedEntries.append( entries) 
                # another way of doing this
                #wantedEntries = pd.concat([wantedEntries, entries])

        # Removing lines that point to same interaction when querying different types of IDs
        filteredTable = wantedEntries.drop_duplicates(subset='interactionID')
        # another way of doing this
        # filteredTable = wantedEntries.copy().groupby(level=0).last()

        print "read_NPInter_file: Number interactions after grepping:",len(filteredTable)

         
        #===================================================================
        # further filtering on NPInter data
        #===================================================================
        # Note: assuming that moleculeB is always the molecule interacting with the RNA
         
        # filter by interaction molecule type ( moleculeBtype must be "protein")
        filteredTable = filteredTable.loc[filteredTable["moleculeBtype"] == KnownScaffoldValidation.MOLECULEBTYPE]
         
        # species must be "Homo sapiens"
        # Note: the species tag seems to be respective of the RNA (or moleculeA), the moleculeB may still be from other species
        filteredTable = filteredTable.loc[filteredTable["species"] == KnownScaffoldValidation.SPECIES] 

        print "read_NPInter_file: Number interactions after molecule type and species filter:",len(filteredTable)
        
        #===================================================================
        # Retrieve set of interacting proteins
        #===================================================================        
        # Note: checking if protein is in RAINET database, to be sure protein is Human, 
        # and to be coherent/fair between catRAPID and NPInter predictions

        wantedColumns = ["moleculeBdatabase","moleculeBID"] #,"moleculeBname"]

        subset = filteredTable.loc[:,wantedColumns]
        tuples = [ tuple(x) for x in subset.values]

        interactingProts = {}

        for tup in tuples:
            proteinDB = tup[0]
            proteinID = tup[1]
            
            if proteinDB == "UniProt":
                if proteinID not in self.proteinsInRainet:
                    # for example this can be protein that belongs to mouse. The previous species filter was relative to the RNA
                    print "read_NPInter_file: ProteinID not found in RAINET database. Using original proteinID: ", proteinID
            else:
                # If database is different than Uniprot, try find uniprotAC using CrossReferences table                
                # lookup ID in crossreferences table and switch to uniprotAC
                if proteinID in self.xrefDict:
                    proteinID = self.xrefDict[proteinID]
                else:
                    print "read_NPInter_file: ProteinID not found in RAINET database, using external source DB. Using original proteinID", proteinID, proteinDB
                             
            if proteinID not in interactingProts:
                interactingProts[proteinID] = 0
            interactingProts[proteinID] += 1

        print "read_NPInter_file: Total number of interacting proteins:",len(interactingProts)

        return interactingProts
 

    # #
    # Reads list of given proteins, returns the ones in RAINET
    # Does xref search and some text processing to improve ID matching
    def read_manual_list_file(self):

        interactingProts = {}

        with open(self.validatedFile,"r") as f:
            listOfProts = [ line.strip() for line in f]

        print "read_manual_list_file: Total number of proteins in file:",len( listOfProts)


        for ID in listOfProts:
            if ID in self.proteinsInRainet:
                proteinID = ID
            else:
                # First: Try to exclude the isoform notation of a uniprotAC # E.g. P09651-2
                tempID = ID.split("-")[0]
                if tempID in self.proteinsInRainet:
                    proteinID = tempID
                else:
                    # Second: try find protein using CrossReferences table                
                    # lookup ID in crossreferences table and switch to uniprotAC
                    if ID in self.xrefDict:
                        proteinID = self.xrefDict[ID]
                    else:
                        print "read_manual_list_file: ProteinID not found in RAINET database. Using original ID:", ID
                        proteinID = ID
                        # continue

            if proteinID not in interactingProts:
                interactingProts[proteinID] = 0
            interactingProts[proteinID] += 1

        print "read_manual_list_file: Total number of interacting proteins:",len( interactingProts)
        
        return interactingProts

    # #
    # Convert the zscore value to interaction propensity
    def convertToPropensity(self, score):
        # catRAPID FAQ: http://s.tartaglialab.com/static_files/shared/faqs.html
        # "Propensity normalised with mean and standard deviation of all fragments. 
        # In the catRAPID omics module, the Z-score is normalized using mean = 23.25 and standard deviation = 37.90 that were calculated on two reference sets. "
        # Z-score = (x - mean) / std , therefore: x = std * z-score + mean
        return score * KnownScaffoldValidation.CATRAPID_ZS_STD + KnownScaffoldValidation.CATRAPID_ZS_MEAN
    

if __name__ == "__main__":
    
    try:
        # Create Logger instance by using the first log action.
        Logger.get_instance().info( "KnownScaffoldValidation : Starting..." )

        #===============================================================================
        # Get input arguments, initialise class
        #===============================================================================
        parser = argparse.ArgumentParser(description='Script to see if known scaffolding RNAs have enriched binding to their known protein targets when using catRAPID data ') 

        # positional args
        parser.add_argument('catRAPIDFile', metavar='catRAPIDFile', type=str,
                             help='File path of CatRAPID omics/fragments results from the webserver.')
        parser.add_argument('validatedFile', metavar='validatedFile', type=str,
                             help='File path of NPInter file or manual list of uniprotACs, one ID per line. Toggle with --npinter flag. E.g. golden_set_NPInter\[v3.0\].txt')
        parser.add_argument('wantedRNAFile', metavar='wantedRNAFile', type=str,
                             help='File path with list of column names and values to search for on NPInter file. This can be any column-value filter. This will retain any line which matches any of these criteria. E.g. moleculeAID\tNONHSAG008670\nmoleculeAName\tNEAT1')    
        parser.add_argument('rainetDB', metavar='rainetDB', type=str, help='Path to RAINET database to be used.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Folder where to write output files.')
        
        # optional args
        parser.add_argument('--discriminativePowerCutoff', metavar='DiscriminativePowerCutoff', default = 0.75, type=float, help='catRAPID Minimum Disciminative power cutoff')
        parser.add_argument('--zscoreCutoff', metavar='zscoreCutoff', default = 0, type=float, help='catRAPID Zscore cutoff')
        parser.add_argument('--topProportion', metavar='topProportion', default = 1, type=float, help='Use float values from 0 to 1. After applying discriminativePower and interactionStrength filters, retrieve given top percent of entries, sorted by interactionStrength and discriminativePower')
        parser.add_argument('--npinter', metavar='npinter', default = 1, type=int, help='Whether validatedFile is NPInter file or simple list of proteins')
        parser.add_argument('--columnForPlot', metavar='columnForPlot', default = 1, choices=[1,2,3], type=int, help='Which column to use for plotting distribution (e.g. 1 for Zscore, 2 for discriminative power)')
        parser.add_argument('--searchSpaceFile', metavar='searchSpaceFile', default = "", type=str, help='File with list of protein uniprotAC to use as search space. Final output will only contain overlap of catRAPID predictions and the proteins in this file.')
        parser.add_argument('--zscoreToPropensity', metavar='zscoreToPropensity', default = 1, type=int, help='Whether to convert Z-score values from input to Interaction propensity, based on static values from catRAPID FAQ.')

    
        #gets the arguments
        args = parser.parse_args( ) 

        # Initialise class
        run = KnownScaffoldValidation( args.catRAPIDFile, args.validatedFile, args.wantedRNAFile, args.rainetDB, 
                                       args.outputFolder, args.discriminativePowerCutoff, args.zscoreCutoff, args.topProportion,
                                       args.npinter, args.columnForPlot, args.searchSpaceFile, args.zscoreToPropensity )

        #===============================================================================
        # Run analysis / processing
        #===============================================================================
         
        # Start chrono
        Timer.get_instance().start_chrono()
 
        Timer.get_instance().step( "reading catRAPID file..")    
        
        # Read CatRAPID
        catRAPIDInteractingProteins, catRAPIDNonInteractingProteins, scoreDict = run.read_catRAPID_file()        

        # Read NPInter or given list of proteins
        
        if run.npinter:
            experimentallyValidatedProteins = run.read_NPInter_file()
        else:
            experimentallyValidatedProteins = run.read_manual_list_file()        

        #===============================================================================
        # Create plot of score with annotation of experimentallyValidatedProteins
        #===============================================================================
        # Note: that currently scores are being picked up before any score filtering
        
        if run.searchSpaceFile != "":
            with open(run.searchSpaceFile,"r") as f:
                searchSpace = { line.strip() for line in f}
        else:
            searchSpace = set( scoreDict.keys())
            assert len( searchSpace.intersection( set(scoreDict.keys())) ) == len( scoreDict)

        print "Length of search space:", len(searchSpace)
        
        outFile = open(run.outputFolder + "/" + run.catRAPIDFile.split("/")[-1]+"_scores.tsv", "w")

        outFile.write("uniprotac\tcatrapid_score\tin_validated_set\n")        
        for prot in scoreDict:
            # Only write to file if protein in given search space (if option provided)
            if prot in searchSpace:
                score = scoreDict[ prot]
                # convert z-score to interaction propensity (if option provided)
                if run.zscoreToPropensity:
                    score = run.convertToPropensity( score)

                outFile.write( "%s\t%s\t%s\n" % ( prot, score, prot in experimentallyValidatedProteins ) )
        
        outFile.close()
        
        # Run R command for creating image
        command = "Rscript %s %s" % ( KnownScaffoldValidation.DISTRIBUTION_SCRIPT, outFile.name)
        result = SubprocessUtil.run_command( command) #, return_stdout = 1, verbose = 1)


        #===============================================================================
        # Hypergeometric test
        #===============================================================================

        # phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
        # 
        # x, q vector of quantiles representing the number of white balls drawn
        # without replacement from an urn which contains both black and white
        # balls.
        # 
        # m the number of white balls in the urn.
        # 
        # n the number of black balls in the urn.
        # 
        # k the number of balls drawn from the urn.

        catRAPIDInteractingSet = set( catRAPIDInteractingProteins.keys())
        catRAPIDNonInteractingSet = set( catRAPIDNonInteractingProteins.keys())

        # background being all catRAPID predictions
        backgroundSet = catRAPIDInteractingSet.union( catRAPIDNonInteractingSet)

        NPInterInteractingSet = set( experimentallyValidatedProteins.keys()).intersection( backgroundSet )               
        notInNPInter = backgroundSet - NPInterInteractingSet

        print ("Experimentally validated Interacting set")
        for key in NPInterInteractingSet:
            print key

        catRAPIDIntNPInter = catRAPIDInteractingSet.intersection( NPInterInteractingSet)

        # white = success
        # black = non-success

        x = len( catRAPIDIntNPInter) # successful draws (i.e. white balls retrieved when drawn)
        m = len( NPInterInteractingSet) # white balls (i.e. total successes in population)
        n = len( notInNPInter) # black balls (i.e. total non-successes in population)
        k = len( catRAPIDInteractingSet) # total draws (i.e. sample size)

        # Run hypergeometric test externally in R and get result
        command = "Rscript %s %s %s %s %s" % ( KnownScaffoldValidation.HYPERGEOMETRIC_TEST_SCRIPT, x, m, n, k)
    
        result = SubprocessUtil.run_command( command, return_stdout = 1, verbose = 0)
        print ("\n")
        print (args)
        print ("Sampling %.2f%% of the population" % ( k * 100.0 / len(backgroundSet)) )
        print ("x: %i\tm: %i\tn: %i\tk: %i" % (x,m,n,k) )
        print ( "Hypergeometric_test p-value:\t%s\n" % result )         
        
        # write pertinent output to file
        outFile = open(run.outputFolder + "/" + run.catRAPIDFile.split("/")[-1], "w")
        outFile.write(str(args)+"\n")
        outFile.write("Sampling %.2f%% of the population\n" % ( k * 100.0 / len(backgroundSet)) )
        outFile.write( "x: %i\tm: %i\tn: %i\tk: %i\n" % (x,m,n,k) )
        outFile.write( "Hypergeometric_test p-value:\t%s\n" % result )         
        outFile.close()


    # Use RainetException to catch errors
    except RainetException as rainet:
        Logger.get_instance().error( "Error during execution of KnownScaffoldValidation. Aborting :\n" + rainet.to_string())

    # Stop the chrono      
    Timer.get_instance().stop_chrono( "KnownScaffoldValidation : Finished" )


#         ## Playing with pandas
# 
#         # Indexing / slicing data frames
#         
#         df = pd.DataFrame({ 1 : [10, 20, 30], 2: [100, 200, 300]})
#         #     1    2
#         # 0  10  100
#         # 1  20  200
#         # 2  30  300
# 
#         print df.loc[1] # return row 1, type Series
#         print "\n"
#         print df[1] # return column 1, type Series
#         print "\n"
#         print df.loc[1,2] # return row 1, column 2, type Numpy Int
#         print "\n"
#         print df[1][2] # return row 2, column 1, type Numpy Int
#         print "\n"
# 
# 
#         print df.loc[:,2] # get column 2, all rows
# 
#         print df.loc[[0,1],[1,2]] # gets rows 0 and 1, columns 1 and 2, type data frame
#         print df.loc[[0,1],[1,10]] # if I invent a column index it will not crash, but return NaN instead
#         print "\n"
# 
#         print df[1][1:3] # gets column 1, rows 1 to 3. type Series
#         print "\n"
# 
#         ## data frame filtering 
# 
#         print df.loc[1] > 20 # return boolean array for row 1
# 
#         print df[1] > 20 # return boolean array for column 1
# 
#         print df.loc[df[1] > 10,:] # return rows for which the value of column 1 is larger than 10
#         print "\n"
# 
#         print df.loc[:, df.loc[0] > 10] # return columns for which the value of row 0 is larger than 10
# # 
#         # Summary: DF loc: df.loc[row,column] xy coordinates. Can use list of items to retrieve smaller DataFrame: df.loc[[0,1],[1,2]]. Cannot use slices
#         #          DF indexing: df[col][num] yx coordinates. Can use slices: df[1][1:3]

