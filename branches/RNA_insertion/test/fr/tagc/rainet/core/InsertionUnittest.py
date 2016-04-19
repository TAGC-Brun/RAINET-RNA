
import unittest
from sqlalchemy import inspect

from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinDomain import ProteinDomain
from fr.tagc.rainet.core.data.ProteinInteraction import ProteinInteraction
from fr.tagc.rainet.core.data.ProteinIsoform import ProteinIsoform
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.KEGGPathway import KEGGPathway
from fr.tagc.rainet.core.data.ProteinKEGGAnnotation import ProteinKEGGAnnotation
from fr.tagc.rainet.core.data.ReactomePathway import ReactomePathway
from fr.tagc.rainet.core.data.ProteinReactomeAnnotation import ProteinReactomeAnnotation
from fr.tagc.rainet.core.data.GeneOntology import GeneOntology
from fr.tagc.rainet.core.data.ProteinGOAnnotation import ProteinGOAnnotation
from fr.tagc.rainet.core.data.NetworkModule import NetworkModule
from fr.tagc.rainet.core.data.NetworkModuleAnnotation import NetworkModuleAnnotation
from fr.tagc.rainet.core.data.PPINetwork import PPINetwork
from fr.tagc.rainet.core.data.PPINetworkInteraction import PPINetworkInteraction
from fr.tagc.rainet.core.data.PartitionAnalysis import PartitionAnalysis
from fr.tagc.rainet.core.data.Gene import Gene
from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.RNACrossReference import RNACrossReference
from fr.tagc.rainet.core.data.MRNA import MRNA
from fr.tagc.rainet.core.data.LncRNA import LncRNA
from fr.tagc.rainet.core.data.OtherRNA import OtherRNA
from fr.tagc.rainet.core.data.ProteinRNAInteractionCatRAPID import ProteinRNAInteractionCatRAPID
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.RNATissueExpression import RNATissueExpression
from fr.tagc.rainet.core.data.Tissue import Tissue

from UnittestConstants import *

# #
# Unittesting a database produced and populated by Rainet InsertionStrategy. 
#
# Note: this is not a test on the Insertion methods themselves but on the output database instead.
# Before running these unittests one needs a Rainet database, populated with test data as in following command:
# Rainet.py Insertion -s human -d /home/diogo/Documents/RAINET_data/TAGC/rainetDatabase/db_testing/rainet_testing_DB.sqlite -i /home/diogo/workspace/tagc-rainet-RNA/resources/insertion_human_rna_test.ini -f
#
class InsertionUnittest(unittest.TestCase):

    # #    
    # Run before each test
    def setUp(self):
        
        self.DBPath = DB_PATH

        # Build a SQL session to DB
        SQLManager.get_instance().set_DBpath(self.DBPath)
        self.sql_session = SQLManager.get_instance().get_session()
        self.db_engine = SQLManager.get_instance().get_engine()


    # #
    # Test correct existence of tables in the database
    def testTables(self):
         
        inspector = inspect(self.db_engine)
                
        tableNames = inspector.get_table_names()
        
        self.assertTrue(len(tableNames) == TOTAL_NUMBER_TABLES,"asserting that number of tables in database is the expected")
        
        for tableName in tableNames:
            columns = inspector.get_columns(tableName)
            if tableName == "RNA":
                self.assertTrue(len(columns) == TOTAL_COLUMNS_IN_RNA_TABLE,"asserting that number of columns in a given table is the expected")
                columnNames = [column['name'] for column in columns]
                self.assertTrue("type" in columnNames,"'type' column represents correspondence to RNA subtype ")


    # #
    # Unittest for RNA SQL table and python objects
    def testRNA(self):
 
        response = self.sql_session.query( RNA ).filter(RNA.transcriptID == EXAMPLE_MRNA).first()
        
        self.assertTrue(isinstance(response, MRNA),"check if the mRNA is instance of MRNA table/class")
        self.assertTrue(isinstance(response, RNA),"check if the mRNA is also instance of RNA table/class")

        definedVars = [var for var in vars(response) if not var.startswith("_")]
        
        self.assertTrue(len(definedVars) == TOTAL_COLUMNS_IN_RNA_TABLE,
                        "asserting that number of attributes corresponds to number of attributes in SQL table.")

        #For a given RNA, see if all values are correct. Note that this also checks if items are string, integer etc
        for var in definedVars:        
            print var    
            value = eval("response."+var)
            self.assertTrue(value == EXAMPLE_MRNA_TABLE_COLUMNS[var], "asserting number of columns in RNA table is correct")
    
            if var == "geneID":
                gene = self.sql_session.query( Gene ).filter(Gene.geneID == value).first()
                self.assertTrue(gene != None,"assert that gene of given RNA is present in Gene table.")
        
        response = self.sql_session.query( MRNA ).filter(MRNA.transcriptID == EXAMPLE_MRNA_PEPTIDE[0]).first()
        self.assertTrue(response.proteinID == EXAMPLE_MRNA_PEPTIDE[1], "asserting insertion of proteinID on the MRNA table")


    # #
    # Unittest for RNACrossReference SQL table and python objects
    def testRNAXRefs(self):

        response = self.sql_session.query( RNACrossReference ).filter(RNACrossReference.transcriptID == EXAMPLE_MRNA).all()
        
        self.assertTrue(len(response) == len(EXAMPLE_MRNA_XREFS), "asserting that number of returned items is correct")
        
        for line in response:
            definedVars = [var for var in vars(line) if not var.startswith("_")]
            
            self.assertTrue(len(definedVars) == TOTAL_COLUMNS_IN_RNAXREF_TABLE,
                            "asserting that number of attributes corresponds to number of attributes in SQL table.")

            tag = "%s|%s" % (line.sourceDB, line.crossReferenceID)
            self.assertTrue(tag in EXAMPLE_MRNA_XREFS, "asserting existence of each cross reference")
                
    # #
    # Unittest for ProteinRNAInteractionCatRAPID SQL table and python objects
    def testPRICatRAPID(self):

        response = self.sql_session.query( ProteinRNAInteractionCatRAPID ).filter(ProteinRNAInteractionCatRAPID.transcriptID == EXAMPLE_PRI_RNA).all()
        
        for line in response:
            
            self.assertTrue(line.peptideID in EXAMPLE_PRI_RNA_INTERACTIONS, "asserting that correct interaction exists")
            
            self.assertTrue(line.interactionScore == EXAMPLE_PRI_RNA_INTERACTIONS[line.peptideID], "asserting that interaction score is correct")

            if line.peptideID in EXAMPLE_PRI_RNA_PROTEIN_UNIPROT:
                self.assertTrue(line.proteinID == EXAMPLE_PRI_RNA_PROTEIN_UNIPROT[line.peptideID],
                                 "asserting that peptideID is associated to correct proteinID")
                protein = self.sql_session.query( Protein ).filter(Protein.uniprotAC == line.proteinID).first()
                self.assertTrue(protein != None,"assert that a given protein in protein-RNA interaction table is present in protein table")
 

    # #
    # Unittest for RNATissueExpression and Tissue SQL table and python objects
    def testRNATissueExpression(self):

        response = self.sql_session.query( RNATissueExpression ).filter(RNATissueExpression.transcriptID == EXAMPLE_MRNA).all()
        tissues = self.sql_session.query( Tissue.tissueName ).all()
        
        tissues = [tiss[0] for tiss in tissues]
        
        self.assertTrue(len(tissues) == NUMBER_TISSUES, "asserting that correct number of tissues exists")
        
        for line in response:
            self.assertTrue( (line.tissueName) in tissues, "asserting that tissues are connected")

        # should add more tests after settling on exact expression data

 
    # #    
    # Run after each test    
    def tearDown(self):
        pass
    

