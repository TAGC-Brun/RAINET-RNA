
import re

from sqlalchemy import or_,and_
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey,PrimaryKeyConstraint
from sqlalchemy.orm import relationship

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger

from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.pattern.PatternUtil import PatternUtil
from fr.tagc.rainet.core.util.sql.SQLUtil import SQLUtil

from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference
from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.Protein import Protein

from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException
from sqlalchemy.orm.base import instance_dict
from fr.tagc.rainet.core.util.time.Timer import Timer
import sys


# #
# This class describes a Protein-RNA interactions predicted by CatRAPID software.
#
class ProteinRNAInteractionCatRAPID( Base ):
    __tablename__ = 'ProteinRNAInteractionCatRAPID'


    # The RNA transcript ID, Ensembl ENST
    transcriptID = Column( String, ForeignKey( 'RNA.transcriptID'), primary_key=True)
    # The Protein Uniprot_ac ID
    proteinID = Column( String, ForeignKey( 'Protein.uniprotAC'), primary_key=True)
    # The ENSEMBL ENSP peptide ID
    peptideID = Column( String, primary_key=True )
    # The interaction score / interaction propensity from catRAPID prediction
    interactionScore = Column( Float )
    # The discriminative power from catRAPID prediction #note: this value is not yet present in testing file
    discriminativePower = Column( Float )
    
    proteins = relationship( "Protein" )

    # #
    # The ProteinRNAInteractionCatRAPID constructor.
    #
    # @param interactors: the interacting protein-RNA pair
    # @param interaction_score: the interaction score
    #
    def __init__( self, interactors, interaction_score):

        from fr.tagc.rainet.core.util.data.DataManager import DataManager

        sql_session = SQLManager.get_instance().get_session()
        dt_manager = DataManager.get_instance()

        # Initialize data items to store missing interactions
        if DataConstants.PROTEIN_RNA_INTERACTION_CATRAPID_MISSING_RNA_KW not in dt_manager.data:
            dt_manager.store_data(DataConstants.PROTEIN_RNA_INTERACTION_CATRAPID_MISSING_RNA_KW,[])
        if DataConstants.PROTEIN_RNA_INTERACTION_CATRAPID_MISSING_PEP_KW not in dt_manager.data:
            dt_manager.store_data(DataConstants.PROTEIN_RNA_INTERACTION_CATRAPID_MISSING_PEP_KW,[])
            
        #=======================================================================
        # Parse interactors
        #
        # note: in the initial catRAPID input, both interactors are in the same TSV column. E.g.  ENSP00000320917_ENST00000389594
        # Peptide is always on left side, RNA in the right side. IDs should not contain "_".
        # Assumption that there only one interaction between each peptide-RNA pair
        #=======================================================================

        spl = interactors.split("_")
        if len(spl) == 2:
            peptide_id = spl[0]
            transcript_id = spl[1]
        else:
            raise RainetException( "ProteinRNAInteractionCatRAPID.__init__ : The interactor string could not be parsed: " + str( interactors ))

        #=======================================================================
        # Fill variables
        #=======================================================================

        try:
            self.interactionScore = float( interaction_score)
        except ValueError as ve:
            raise RainetException( "ProteinRNAInteractionCatRAPID.__init__ : The value of interaction score is not a float: " + str( interaction_score ), ve )

        #Store peptide_id so that we can retrace the isoform used in catRAPID
        self.peptideID = peptide_id

        #=======================================================================
        # Query Protein uniprot_ac using cross references
        # From Cross ref ID, search primary ID using cross reference DataManager item created during insertion
        # Note: assuming that if uniprot_ac exists in protein cross reference table it will also exist in protein table.
        #=======================================================================
 
        proteinXrefs = dt_manager.get_data( DataConstants.PROTEIN_ENSP_XREF_KW)

        if peptide_id in proteinXrefs:
            self.proteinID = proteinXrefs[ peptide_id][0]
        else:
#            raise RainetException( "ProteinRNAInteractionCatRAPID.init : No Protein object while using cross references for peptide_id = " + peptide_id)
            Logger.get_instance().warning( "\nProteinRNAInteractionCatRAPID.init : Peptide ID not found, will skip interaction:\t" + str( peptide_id) )
            # Store missing peptide ID in a list
            dt_manager.data[ DataConstants.PROTEIN_RNA_INTERACTION_CATRAPID_MISSING_PEP_KW].append( peptide_id)
            raise NotRequiredInstantiationException( "ProteinRNAInteractionCatRAPID.init : No Protein found, instance will not be created.")
 
        #=======================================================================
        # Query RNA object
        # See if RNA with given transcript_id exists in database
        #=======================================================================

        RNA_list = dt_manager.get_data( DataConstants.RNA_ALL_KW)
 
        if transcript_id in RNA_list:
            self.transcriptID = transcript_id
        else:
#            raise RainetException( "ProteinRNAInteractionCatRAPID.init : No RNA object found for transcript_id = " + transcript_id)
            Logger.get_instance().warning( "\nProteinRNAInteractionCatRAPID.init : RNA ID not found, will skip interaction:\t" + str( transcript_id) )
            # Store missing RNA ID in a list
            dt_manager.data[ DataConstants.PROTEIN_RNA_INTERACTION_CATRAPID_MISSING_RNA_KW].append( transcript_id)
            raise NotRequiredInstantiationException( "ProteinRNAInteractionCatRAPID.init: No RNA found, instance will not be created." )


    ##
    # Add the object to SQLAlchemy session if it is linked to a protein and RNA
    def add_to_session(self):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self)
