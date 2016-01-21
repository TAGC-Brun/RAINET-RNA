
from sqlalchemy import Column, String, Integer, Boolean
from sqlalchemy.orm import relationship

from sets import Set

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger


from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.data.GeneSymbol import GeneSymbol
from fr.tagc.rainet.core.data.SynonymGeneSymbol import SynonymGeneSymbol
from fr.tagc.rainet.core.data.ProteinDomain import ProteinDomain
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.ProteinIsoform import ProteinIsoform
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from sqlalchemy.orm.base import instance_dict



# #
# This class describe a protein obtained from query on the Uniprot DB
#
# this class has a one-to-many relationship with the class SynonymGeneSymbol
#
class Protein( Base ):
    __tablename__ = 'Protein'

    # The Uniprot accession number of the protein
    uniprotAC = Column( String, primary_key = True )
    # The Uniprot ID of the protein
    uniprotID = Column( String )
    # The Uniprot ID of the protein
    uniprotName = Column( String )
    # The organism of the protein
    organism = Column( String )
    # The length of protein
    proteinLength = Column( Integer )
    # The ???????????
    idFragment = Column( Boolean )
    # The list of Uniprot Gene Symbol
    uniprotGeneSymbols = relationship( 'GeneSymbol', backref = 'protein' )
    # The list of protein domains
    domains = relationship( 'ProteinDomain', backref = 'protein' )
    # The list of protein domains
    crossReferences = relationship( 'ProteinCrossReference', backref = 'protein' )
    # The list of protein isoforms
    uniprotIsoforms = relationship( 'ProteinIsoform', backref = 'protein' )
    # The list of ProteinInteractions
    proteinInteractions = relationship( 'ProteinInteraction')
    
    # #
    # The Protein constructor
    # 
    # @param prot_ac : string - The Uniprot accession number of the protein
    # @param prot_id : string - The Uniprot ID of the protein
    # @param prot_name : string - The Uniprot ID of the protein
    # @param prot_symbol : string - The (list of) Uniprot gene symbols of the protein
    # @param prot_synonym : string - The (list of) Uniprot gene symbol synonyms of the protein
    # @param prot_organism : string - The organism of the protein
    # @param prot_length : string - The length of protein
    # @param prot_fragment : boolean - 
    def __init__( self, prot_ac, prot_id, prot_name, prot_symbol, prot_synonym, prot_organism, prot_length, prot_fragment, prot_domain_pfam, prot_domain_smart ):
        
        #=======================================================================
        # Fill the main protein variables
        #=======================================================================
        self.uniprotAC = prot_ac
        self.uniprotID = prot_id
        self.uniprotName = prot_name
        self.organism = prot_organism
        try:
            self.proteinLength = int( prot_length )
        except ValueError as ve:
            raise RainetException( "Protein.__init__ : The value of protein length is not an integer: " + str( prot_length ), ve )
            
        if prot_fragment == None or prot_fragment == '':
            self.idFragment = False
        else:
            self.idFragment = True
        
        #=======================================================================
        # Build the GeneSymbol objects related to the Protein and their 
        # corresponding SynonymeGeneSymbol objects
        #
        # Note that gene symbols are separated by ";" in the prot_symbol term
        # and the synonyms are separated by ";" in prot_synonyms trem 
        # when they are from different gene symbol and by " " when they are 
        # from the same gene symbol.
        # for instance:
        # prot_symbol = "USP17L24; USP17L25; USP17L26; USP17L27; USP17L28; USP17L29; USP17L30"
        # prot_synonym = "USP17 USP17H USP17I USP17J USP17K USP17L USP17M; ; ; ; ; ;"
        # means "USP17L24" have 7 synonyms : USP17 USP17H USP17I USP17J USP17K USP17L USP17M
        # whereas the other symbols have no synonyms
        #=======================================================================
        self.uniprotGeneSymbols = []
        if prot_symbol != None and prot_symbol != "":
            # Get the list of symbols and the list of corresponding synonyms
            symbol_list = prot_symbol.split( ";" )
            synonym_tag_list = prot_synonym.split( ";" )
            for current_index in range( 0, len( symbol_list ) ):
                current_symbol = symbol_list[ current_index]
                if current_symbol != None and current_symbol != "":
                    # Build the new GeneSymbol and add it to the list of protein symbols 
                    new_symbol = GeneSymbol( uniprotGeneSymbol = current_symbol.strip() )
                    self.uniprotGeneSymbols.append( new_symbol )
                    # Build the SynonymGeneSymbols associated to the GeneSymbol if any
                    if len( synonym_tag_list ) >= ( current_index + 1 ):
                        current_synonym_tag = synonym_tag_list[current_index]
                        if current_synonym_tag != None and current_synonym_tag != "":
                            current_synonym_list = set( current_synonym_tag.split( " " ))
                            for current_synonym in current_synonym_list:
                                if current_synonym != None and current_synonym != "":
                                    new_synonym = SynonymGeneSymbol( uniprotSynonymGeneSymbol = current_synonym.strip())
                                    new_symbol.add_synonym( new_synonym)
                        
        
        #=======================================================================
        # Build the ProteinDomain objects related to the Protein
        #=======================================================================
        # -- Build the domains from PFAM
        if prot_domain_pfam != None and prot_domain_pfam != "":
            domain_list = prot_domain_pfam.split( ";" )
            for domain in domain_list:
                if domain != None and domain != '':
                    try:
                        new_domain = ProteinDomain( domain, "", DataConstants.PROTEIN_DOMAIN_DBSOURCE_PFAM )
                        self.domains.append( new_domain )
                    except NotRequiredInstantiationException:
                        pass
        # -- Build the domains from SMART
        if prot_domain_smart != None and prot_domain_smart != "":
            domain_list = prot_domain_smart.split( ";" )
            for domain in domain_list:
                if domain != None and domain != '':
                    try:
                        new_domain = ProteinDomain( domain, "", DataConstants.PROTEIN_DOMAIN_DBSOURCE_SMART )
                        self.domains.append( new_domain )
                    except NotRequiredInstantiationException:
                        pass

        
    # #
    # Add the object to SQLAlchemy session
    def add_to_session( self ):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self )

    # #
    # Add a ProteinCrossReference to the protein cross reference list
    #
    # @param cross_reference : string - A ProteinCrossReference instance
    def add_cross_reference( self, cross_reference ):
        
        if cross_reference != None:
            self.crossReferences.append( cross_reference )
        
    # #
    # Add a ProteinIsoform to the protein isoform list
    #
    # @param isoform : string - A ProteinIsoform instance
    def add_isoform( self, isoform ):
        
        if isoform != None:
            self.uniprotIsoforms.append( isoform )
