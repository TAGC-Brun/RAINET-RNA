from sqlalchemy.sql.schema import Column
from sqlalchemy.sql.sqltypes import String

from fr.tagc.rainet.core.util.sql.Base import Base

# #
# This class represent the status of the DB table, indicating if the table was correctly filled,
# filled with warning or with error or not filled
# 
class TableStatus( Base ):
    __tablename__ = "TableStatus"
    
    # The name of the table
    tableName = Column( String, primary_key = True )
    # The Table status
    tableStatus = Column( String )
    # The table source (file used to fill the table)
    tableSource = Column( String )
    
    
    def __init__(self, name, status, source):
        
        self.tableName = name
        self.tableStatus = status
        self.tableSource = source
        
    
        
        
    
