
from sqlalchemy import Column, String, Integer, Boolean, Float, ForeignKey

from fr.tagc.rainet.core.data.RNA import RNA
from fr.tagc.rainet.core.data.Protein import Protein
from fr.tagc.rainet.core.data.ProteinCrossReference import ProteinCrossReference

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.log.Logger import Logger


class MRNA( RNA ):
    
    __tablename__ = 'MRNA'
    
    transcriptID = Column( String, ForeignKey('RNA.transcriptID'), primary_key=True)

    # The Protein Uniprot_ac ID
    proteinID = Column( String, ForeignKey( 'Protein.uniprotAC'))


    __mapper_args__ = { 'polymorphic_identity':'MRNA'}



    def __init__( self, peptide_id ):
        
        #=======================================================================
        # Query Protein uniprot_ac using cross references
        # From Cross ref ID, search primary ID using cross reference DataManager item created during insertion
        # Note: assuming that if uniprot_ac exists in protein cross reference table it will also exist in protein table.
        #=======================================================================
 
        from fr.tagc.rainet.core.util.data.DataManager import DataManager

        proteinXrefs = DataManager.get_instance().get_data(DataConstants.PROTEIN_ENSP_XREF_KW)
        
        if peptide_id in proteinXrefs:
            self.proteinID = proteinXrefs[peptide_id][0]
        else:
            Logger.get_instance().debug( "\n MRNA.init : Peptide ID not found:\t" + str(peptide_id) )

