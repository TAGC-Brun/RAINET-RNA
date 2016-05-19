
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
    # The interaction score / interaction propensity from catRAPID prediction
    interactionScore = Column( Float )
    
    proteins = relationship( "Protein" )

    # #
    # The ProteinRNAInteractionCatRAPID constructor.
    #
    # @param interactors: the interacting protein-RNA pair
    # @param interaction_score: the interaction score
    #
    def __init__( self, interactors, interaction_score):

        from fr.tagc.rainet.core.util.data.DataManager import DataManager

        dt_manager = DataManager.get_instance()
    
        #=======================================================================
        # Parse interactors
        #
        # Example
        # sp|Q96DC8|ECHD3_HUMAN ENST00000579524   -12.33  0.10    0.00
        # sp|P10645|CMGA_HUMAN ENST00000516610    10.66   0.32    0.00
        # protein and rna separated by " ", other values separated by "\t"
        # 
        # Protein is always on left side, RNA in the right side.
        # Assumption that there only one interaction between each Protein-RNA pair
        #=======================================================================

        spl = interactors.split(" ")
        if len(spl) == 2:
            protein_id = spl[0].split( "|")[1]
            transcript_id = spl[1].split( "\t")[0]
        else:
            raise RainetException( "ProteinRNAInteractionCatRAPID.__init__ : The interactor string could not be parsed: " + str( interactors ))

        #=======================================================================
        # Fill variables
        #=======================================================================

        try:
            self.interactionScore = float( interaction_score)
        except ValueError as ve:
            raise RainetException( "ProteinRNAInteractionCatRAPID.__init__ : The value of interaction score is not a float: " + str( interaction_score ), ve )

        #=======================================================================
        # Query Protein object
        # See if Protein with given protein_id exists in database
        #=======================================================================
 
        protein_list = dt_manager.get_data( DataConstants.PROT_ALL_KW)

        if protein_id in protein_list:
            self.proteinID = protein_id
        else:
#            raise RainetException( "ProteinRNAInteractionCatRAPID.init : No Protein object while using cross references for protein_id = " + protein_id)
            Logger.get_instance().warning( "\nProteinRNAInteractionCatRAPID.init : Protein ID not found, will skip interaction:\t" + str( protein_id) )
            # Store missing Protein ID in a list
            dt_manager.data[ DataConstants.PROTEIN_RNA_INTERACTION_CATRAPID_MISSING_PROT_KW].append( protein_id)
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
