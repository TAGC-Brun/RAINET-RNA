
from sqlalchemy import Column, String, Boolean, ForeignKey,PrimaryKeyConstraint

from fr.tagc.rainet.core.util.sql.Base import Base

from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException

##
# This class describe an isoform associated to a protein
# this class has a many-to-one relationship with the class Protein
#
class ProteinIsoform( Base):
    __tablename__ = 'ProteinIsoform'
    
    # The base Protein
    protein_id = Column( String, ForeignKey('Protein.uniprotAC'))
    # The synonym symbol
    isoformUniprotAC = Column( String)
    # The indication if isoform is a canonical one or an alternative one
    isAlternative = Column( Boolean)
    # Define the Composite PrimaryKey
    __table_args__ = (
        PrimaryKeyConstraint('protein_id', 'isoformUniprotAC'),
    )
    
    ##
    # The constructor of the class
    #
    # @param iso_ac : string - the isoform UniprotAC (for instance P26510-2)
    # @param is_alternative : boolean - False if the isoform is a canonical form, True if not
    def __init__(self, iso_ac, protein_ac, is_alternative):
        
        # Assign the basic parameters
        self.isoformUniprotAC = iso_ac
        self.isAlternative = bool( is_alternative)
        
        # Found the uniprotAC in the isoformUniprotAC as the string before the "-" character
        # for instance : isoformUniprotAC=P65321-2 => uniprotAC=P65321
        if iso_ac != None and iso_ac != "":
            if protein_ac == None:
                raise RainetException("ProteinIsoform.__init__: abnormal isoform definition : Isoform AC = " + iso_ac + " but original protein uniprotAC is not defined.")
            # Get the Protein corresponding to the given unirptoAC
            sql_session = SQLManager.get_instance().get_session()
            from fr.tagc.rainet.core.data.Protein import Protein
            protein_list = sql_session.query( Protein).filter( Protein.uniprotAC == protein_ac).all()
            # If the protein list is correct (one result), build the relationship between the Protein and the ProteinIsoform
            if protein_list != None and len( protein_list) != 0:
                if len( protein_list) == 1:
                    protein = protein_list[0]
                    if protein != None:
                        protein.add_isoform( self)
                    else:
                        raise RainetException("ProteinIsoform.__init__: Found Proteins with uniprotAC= " + protein_ac + " is None.")
                else:
                    raise RainetException("ProteinIsoform.__init__: Several Proteins found with uniprotAC= " + protein_ac + " : " + str( len( protein_list)) + " Proteins")
            else:
                raise NotRequiredInstantiationException("ProteinIsoform.__init__: No Protein found with uniprotAC= " + protein_ac)
        else:
            raise RainetException("ProteinIsoform.__init__: No isoformProteinAC provided: ")

    # #
    # Add the object to SQLAlchemy session
    def add_to_session( self ):
    
        sql_session = SQLManager.get_instance().get_session()
        sql_session.add( self )