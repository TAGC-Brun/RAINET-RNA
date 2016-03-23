
import sys
import os
import argparse
import glob
import numpy as np
import random 
import pandas as pd

from fr.tagc.rainet.core.util.file.FileUtils import FileUtils
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger
from fr.tagc.rainet.core.util.time.Timer import Timer
from fr.tagc.rainet.core.util.subprocess.SubprocessUtil import SubprocessUtil
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinInteraction import ProteinInteraction


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
# - Read supposed interacting proteins from NPInter database file.
# - Match those proteins to RAINET database using external Refs
# - Read catRAPID omics/fragments file, apply interaction cutoff, retrieve list of proteins (uniprotAC/ID) that are in RAINET database, 
# - Perform enrichment analysis 
#===============================================================================

#===============================================================================
# Processing notes:
# 
# Assuming that NPInter file has RNA as moleculeA and partner as moleculeB
#
#===============================================================================


class KnownScaffoldValidation( object ):

    SPECIES = "Homo sapiens"
    MOLECULEBTYPE = "protein"

    def __init__(self, catRAPIDFile, npinterFile, wantedRNAFile, rainetDB, outputFolder, discriminativePowerCutoff, interactionStrengthCutoff):

        self.catRAPIDFile = catRAPIDFile
        self.npinterFile = npinterFile
        self.wantedRNAFile = wantedRNAFile
        self.rainetDB = rainetDB
        self.outputFolder = outputFolder
        self.discriminativePowerCutoff = discriminativePowerCutoff
        self.interactionStrengthCutoff = interactionStrengthCutoff

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.rainetDB)
        self.sql_session = SQLManager.get_instance().get_session()

        self.proteinsInRainet = self.proteins_in_rainet()


    # #
    # Access RAINET database and retrieve all uniprotACs included in the database
    def proteins_in_rainet(self):
        
        query = self.sql_session.query( Protein.uniprotAC ).all()
        
        proteinsInRainet = { str(prot[0]) for prot in query}

        return proteinsInRainet

    # #
    # Read catRAPID file, apply wanted interaction score filters and retrieve set of interacting proteins
    def read_catRAPID_file(self):

        ## Protein ID   Transcript ID   Z-score   Discriminative Power   Interaction Strength   Motif Reference   Domain Presence   Domain Score   -   Star Rating Score
        # Note: catRAPID file header does not match content and is not tab-delimited

        # read file skipping header line 
        table = pd.read_table( self.catRAPIDFile, header = None, sep = "\t", skip_blank_lines = True, skiprows = 1)

        print "Starting number of interactions:",len(table)               
#         print table.head(5)
        
        # #
        # Filter for discriminative power
        table = table.loc[table.loc[:,2] > self.discriminativePowerCutoff,:]

        print "Number interactions after discriminative power filter:",len(table)

        # #
        # Filter for interaction strength
        table = table.loc[table.loc[:,3] > self.interactionStrengthCutoff,:]

        print "Number interactions after interaction strength filter:",len(table)

        # # 
        # get set of interacting proteins
        
        # Note: the transcript may have been fragmented, so there will be several hits for same protein - transcript pair.
        setOfProteins = {}

        # get item storing transcript-protein pair
        interactions = table.loc[:,0]
        
        for inter in interactions:
            # split different protein IDs and transcript
            spl = inter.split("|")
            uniprotAC = spl[1]
            
                        
            if uniprotAC in self.proteinsInRainet:
                if uniprotAC not in setOfProteins:
                    setOfProteins[uniprotAC] = 0
                    
                setOfProteins[uniprotAC] += 1
            else:
                # this can be as proteinAC was deprecated etc
                print "read_catRAPID_file: ProteinID not found in RAINET database", uniprotAC

        print "read_catRAPID_file: Total number of interacting proteins:",len(setOfProteins)
        
        return setOfProteins


    # #
    # Read NPInter file and retrieve list of proteins interacting with wanted RNA
    def read_NPInter_file(self):
        
        # Note: NPInter uses NONCODE database for transcript IDs, but in fact uses Gene IDs, not the transcript IDs
        
        
        # #
        # read wantedRNAFile,  file containing IDs to look for on NPInter file
        wantedRNATable = pd.read_table( self.wantedRNAFile, header = None, sep = "\t", skip_blank_lines = True)

        
        # #
        # read NPInter file using header line 
        table = pd.read_table( self.npinterFile, header = 0, sep = "\t", skip_blank_lines = True)
 
        # #
        # get lines containing wanted RNA in wanted column
        
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

        filteredTable = wantedEntries.drop_duplicates(subset='interactionID')
        # another way of doing this
        # filteredTable = wantedEntries.copy().groupby(level=0).last()

         
        # #
        # further filtering
         
        # filter by interaction molecule type ( moleculeBtype must be "protein")
        filteredTable = filteredTable.loc[filteredTable["moleculeBtype"] == KnownScaffoldValidation.MOLECULEBTYPE]
         
        # species must be "Homo sapiens"
        filteredTable = filteredTable.loc[filteredTable["species"] == KnownScaffoldValidation.SPECIES]
  

        # Note: assuming that moleculeB is always the molecule interacting with the RNA
        wantedColumns = ["moleculeBdatabase","moleculeBID"] #,"moleculeBname"]

        subset = filteredTable.loc[:,wantedColumns]
        tuples = [ tuple(x) for x in subset.values]
            #loc[:,wantedColumns]

        setOfProteins = {}

        for tup in tuples:
            proteinDB = tup[0]
            proteinID = tup[1]
            
            if proteinDB == "UniProt":
                query = self.sql_session.query( Protein ).filter( Protein.uniprotAC == proteinID).all()
            
                if len(query) > 0:
                    if proteinID not in setOfProteins:
                        setOfProteins[proteinID] = 0
                    setOfProteins[proteinID] += 1                    
                else:
                    # for example this can be protein that belongs to mouse. The previous species filter was relative to the RNA
                    print "read_NPInter_file: ProteinID not found in RAINET database", proteinID
            else:
                # TODO: add check of cross references
                print (proteinDB)
        
        print "read_NPInter_file: Total number of interacting proteins:",len(setOfProteins)
        return setOfProteins


        

if __name__ == "__main__":
    
    try:

        #===============================================================================
        # Get input arguments, initialise class
        #===============================================================================
        parser = argparse.ArgumentParser(description='Script to see if known scaffolding RNAs have enriched binding to their known protein targets when using catRAPID data ') 

        # positional args
        parser.add_argument('catRAPIDFile', metavar='catRAPIDFile', type=str,
                             help='CatRAPID omics/fragments results from the webserver.')
        parser.add_argument('npinterFile', metavar='npinterFile', type=str,
                             help='File from NPInter. E.g. golden_set_NPInter\[v3.0\].txt')
        parser.add_argument('wantedRNAFile', metavar='wantedRNAFile', type=str,
                             help='File with list of column names and values to search for on NPInter file. This can be any column-value filter. This will retain any line which matches any of these criteria. E.g. moleculeAID\tNONHSAG008670\nmoleculeAName\tNEAT1')    
        parser.add_argument('rainetDB', metavar='rainetDB', type=str, help='Path to RAINET database to be used.')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Folder where to write output files.')
        
        # optional args
        parser.add_argument('--discriminativePowerCutoff', metavar='DiscriminativePowerCutoff', default = 0.5, type=float, help='catRAPID Minimum Disciminative power cutoff')
        parser.add_argument('--interactionStrengthCutoff', metavar='InteractionStrengthCutoff', default = 0.5, type=float, help='catRAPID Interaction Strength cutoff')
        
        #display help when misusage
        if len(sys.argv) < 5: 
            parser.print_help()
    
        #gets the arguments
        args = parser.parse_args( ) 
        
        # Initialise class
        run = KnownScaffoldValidation( args.catRAPIDFile, args.npinterFile, args.wantedRNAFile, args.rainetDB, 
                                       args.outputFolder, args.discriminativePowerCutoff,
                                        args.interactionStrengthCutoff )

        #===============================================================================
        # Run analysis / processing
        #===============================================================================
         
        # Create Logger instance by using the first log action.
        Logger.get_instance().info( "KnownScaffoldValidation : Starting..." )
 
        # Start chrono
        Timer.get_instance().start_chrono()
 
        Timer.get_instance().step( "reading catRAPID file..")    
        
        catRAPIDInteractingProteins = run.read_catRAPID_file()
        
        NPInterInteractingProteins = run.read_NPInter_file()
        
        catRAPIDSet = set(catRAPIDInteractingProteins.keys())
        
        NPInterSet = set(NPInterInteractingProteins.keys())
        
        print catRAPIDSet.intersection(NPInterSet)
        
        ### CONTINUE HERE!! NEED TO MAKE STATISTICAL TEST
        

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
