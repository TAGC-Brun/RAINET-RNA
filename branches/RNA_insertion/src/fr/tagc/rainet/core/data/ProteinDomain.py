
from sqlalchemy import Column, String, Integer, ForeignKey, PrimaryKeyConstraint
from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql.SQLManager import SQLManager
from fr.tagc.rainet.core.util.exception.NotRequiredInstantiationException import NotRequiredInstantiationException

# #
# This class describes a protein domain
#
class ProteinDomain( Base ):
    __tablename__ = 'ProteinDomain'
    
    # The base Protein
    protein_id = Column( String, ForeignKey( 'Protein.uniprotAC' ) )
    # The domain ID
    domainID = Column( String )
    # The database the information was obtained from
    sourceDB = Column( String )
    # The domain name
    domainName = Column( String )
    # Define the Composite PrimaryKey
    __table_args__ = (
        PrimaryKeyConstraint('protein_id', 'domainID', "sourceDB"),
    )
    
    #
    # The constructor of the class
    #
    # @param domain_id : string - The ID of the domain
    # @param domain_name : string - The name of the domain
    # @param db_source : string - The database name from which the domain was retrieved
    def __init__(self, domain_id, domain_name, db_source):
        
        # Get a SQLalchemy session
        sql_session = SQLManager.get_instance().get_session()
        
        # Retrieve the list of ProteinDomain already in database with the same domainID and DBSource as provided
        domain_list = sql_session.query( ProteinDomain).filter( ProteinDomain.domainID == domain_id, ProteinDomain.sourceDB == db_source).all()

        # If some ProteinDomain already exists in database, update them with the value of domainName
        # and raise a NotRequiredInstantiationException to indicate no new object is required
        if domain_list != None and len( domain_list) != 0:
            for protein_domain in domain_list:
                protein_domain.domainID = domain_id
                protein_domain.domainName = domain_name
                protein_domain.sourceDB = db_source
                sql_session.add( protein_domain)
            raise NotRequiredInstantiationException( "ProteinDomain.init : ProteinDomain already exists with" + domain_id + " and sourceDB=" + db_source)
        
        # At this point, it means there is no corresponding ProteinDomain in the database,
        # so fill the new object variables with the correct values 
        self.domainID = domain_id
        self.domainName = domain_name
        self.sourceDB = db_source
    
    ##
    # Add the object to SQLAlchemy session if it is linked to a protein
    def add_to_session(self):
    
        if self.protein_id != None and self.protein_id != '':
            sql_session = SQLManager.get_instance().get_session()
            sql_session.add( self)
    
