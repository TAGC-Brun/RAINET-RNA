# -*- coding: utf-8 -*-

from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from sqlalchemy import exc
from posix import remove

from fr.tagc.rainet.core.util.log.Logger import Logger

from fr.tagc.rainet.core.util.sql.Base import Base
from fr.tagc.rainet.core.util.sql import SQLConstants

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.data import DataConstants
from fr.tagc.rainet.core.util.option import OptionConstants
from fr.tagc.rainet.core.util.option.OptionManager import OptionManager
from fr.tagc.rainet.core.data.DBParameter import DBParameter

## This class is a singleton aiming to manage SQL connection to the database
# The singleton is able to manage the creation of a SQLAlchemy Session to the database
# it has been initiated with. The singleton keep the same session open until it is
# asked to close it.
class SQLManager( object ) :

    __instance = None

    def __init__( self ):
        self.DBPath = None
        self.session = None


    # # 
    # Returns the current SQLAlchemy session if it exists or create a new one
    # if no previous session was opened
    #
    # @return a new or the current SQLAlchemy session
    def get_session( self ):

        if( self.session == None ):
            # Create engine to dedicated database
            engine = create_engine( SQLConstants.PATH_SQLALCHEMY_SQLITE + self.DBPath )

            # Open the DB session
            session = sessionmaker()
            session.configure( bind = engine, autoflush = True, expire_on_commit = False )

            # Get the session and insert several objects in ENSEMBL table
            self.session = session()

        return self.session

    # #
    # Close the session and insure new session will be created
    # when SQLsession will be requested next time
    #
    # @return None
    def close_session( self ):
        
        self.session.close()
        self.session = None
        return
    
    # #
    # Commit the actual session and close it after
    #
    # @return None
    def commit( self ):
        
        try:
            self.session.commit()
        except exc.SQLAlchemyError as sqle:
            self.rollback_session()
            raise RainetException( "SQLManager.commit : An error occurred while committing the session.", sqle )
        finally:
            self.close_session()

    # #
    # Rollback the actual session
    #
    # @return None
    def rollback_session(self):
        
        try:
            self.session.rollback()
        except exc.SQLAlchemyError as sqle:
            raise RainetException( "SQLManager.rollback : An error occurred while rollbacking the session.", sqle )
    
    # #
    # Set the Database path used by the SQLAlchemu session
    #
    # @param path : string - The path to the database
    #
    # @return None
    def set_DBpath( self, path ):

        self.DBPath = path

    # #
    # Remove the database file at the given path
    # 
    # @param path : string - The path to the database file
    def remove_database_file( self, path ):
        
        try:
            remove( path )
            Logger.get_instance().info( 'File : "' + str( path ) + '" deleted.' )
        except:
            Logger.get_instance().info( 'Can not delete file : ' + str( path ) )
            pass
        

    # #
    # Build the database schema on the provided database path
    # 
    # @param path : string - The path to the database
    # @param keep_file : boolean - Indicates if the previous DB file should be kept or not
    #
    # @return None
    def build_database( self, path, force_override ):
        
        # Get the value of database species specified by the user
        species = OptionManager.get_instance().get_option( OptionConstants.OPTION_SPECIES)
        if species == None or len( species) == 0:
            raise RainetException( "SQLManager.build_database: You must specify a species in your command line. Please check help.")
        
        # Remove the database file if required
        if force_override:
            self.remove_database_file( path )

        # Create engine to dedicated database
        engine = create_engine( SQLConstants.PATH_SQLALCHEMY_SQLITE + path )

        # Open the DB session
        session = sessionmaker()
        session.configure( bind = engine )

        # Look if the Protein table exists. If not, it means DB model was not created
        model_exists = True
        try:
            table_content = session().query( DataConstants.PROTEIN_CLASS).all()
            if table_content == None:
                model_exists = False
        except Exception:
            model_exists = False

        # Create all the required table in DB according to class model
        # if required (override mode or model does not exist in DB)        
        if force_override or not model_exists:
            Base.metadata.create_all( engine )
            sql_session = session()
            SQLManager.check_species(sql_session, species, True)
            sql_session.close()

        # Keep the DB path
        self.DBPath = path

        Logger.get_instance().info( 'Database File created/used : ' + str( path ) )
            

    #
    # Check if the database has been defined using the species specified by the user
    #
    # @return None
    # @throw RainetException when the species defined in the databse does not correspond
    # to the user specified species
    @staticmethod
    def check_species( sql_session, species, create = False):
        
        try:
            db_species = sql_session.query( DBParameter.parameterValue).filter( DBParameter.parameterName == OptionConstants.OPTION_SPECIES).first()
            if db_species == None or len( db_species) == 0:
                if create:
                    species_param = DBParameter( OptionConstants.OPTION_SPECIES, "string", species)
                    sql_session.add( species_param)
                    sql_session.commit()
                else:
                    raise RainetException( "SQLManager.check_species: The database has no species declared and must not be used with that command.")
            else:
                db_species_string = db_species[0]
                if db_species_string.lower() != species.lower():
                    raise RainetException( "SQLManager.check_species: The database species is not the same as the user's specified species: " + db_species_string + " != " + species)
        except exc.SQLAlchemyError as sqle:
            if create:
                sql_session.rollback()
                raise RainetException( "SQLManager.build_databse : An error occurred while committing the species parameter.", sqle )
            else:
                raise RainetException( "SQLManager.build_databse : An error occurred while querying the species parameter.", sqle )
        
        
    # #
    # Returns the singleton instance
    #
    # @return the singleton instance
    @staticmethod
    def get_instance():

        if SQLManager.__instance == None:
            SQLManager.__instance = SQLManager()
        return SQLManager.__instance

    # #
    # Returns the SQLAlchemy engine of the DB
    # 
    # @return the SQLAlchemy engine of the DB
    def get_engine( self):

        if self.DBPath != None:
            engine = create_engine( SQLConstants.PATH_SQLALCHEMY_SQLITE + self.DBPath)
            return engine
        
        return None
        
