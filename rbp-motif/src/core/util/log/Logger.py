# -*- coding: utf-8 -*-

import logging
from logging.handlers import RotatingFileHandler

from core.util import Constants

##
# This class is a singleton used to manage logging
# It contains several method to add log at different level and
# manage both output to file and standard output
class Logger(object):

    __instance = None

    def __init__(self):
        self.logg = self.setLogger()


    ##
    #
    #
    def setLogger( self):
        # Reinitialize log file
        ERROR_FILE = open(Constants.PATH_LOG, 'a')
        ERROR_FILE.write("")
        ERROR_FILE.close()

        # création de l'objet logger qui va nous servir à écrire dans les logs
        logger = logging.getLogger()
        # on met le niveau du logger à DEBUG, comme ça il écrit tout
        logger.setLevel(logging.DEBUG)

        # création d'un formateur qui va ajouter le temps, le niveau
        # de chaque message quand on écrira un message dans le log
        formatter = logging.Formatter \
        ('%(asctime)s :: %(levelname)s :: %(message)s')

        # création d'un handler qui va rediriger une écriture du log vers
        # un fichier en mode 'append', avec 100 backup et une taille max de 50Mo
        self.fileHandler = RotatingFileHandler(Constants.PATH_LOG, 'a', 50000000, 100)

        # Set level on DEBUG
        # créé précédement et on ajoute ce handler au logger
        self.fileHandler.setLevel(logging.DEBUG)
        self.fileHandler.setFormatter(formatter)
        logger.addHandler(self.fileHandler)

        # Création d'un second handler qui va rediriger chaque écriture de log
        # sur la console
        self.streamHandler = logging.StreamHandler()
        self.streamHandler.setLevel(logging.DEBUG)
        logger.addHandler(self.streamHandler)

        return logger

    ##
    # Set the level of verbosity to the given level
    # 
    # @param level : string - The level of verbosity
    def set_level(self, level):
        
        if level == Constants.MODE_DEBUG:
            self.fileHandler.setLevel( logging.DEBUG)
            self.streamHandler.setLevel( logging.DEBUG)
        elif level == Constants.MODE_INFO:
            self.fileHandler.setLevel( logging.INFO)
            self.streamHandler.setLevel( logging.INFO)
        elif level == Constants.MODE_WARNING:
            self.fileHandler.setLevel( logging.WARNING)
            self.streamHandler.setLevel( logging.WARNING)
        elif level == Constants.MODE_CRITICAL:
            self.fileHandler.setLevel( logging.CRITICAL)
            self.streamHandler.setLevel( logging.CRITICAL)
        elif level == Constants.MODE_ERROR:
            self.fileHandler.setLevel( logging.ERROR)
            self.streamHandler.setLevel( logging.ERROR)
        elif level == Constants.MODE_FATAL:
            self.fileHandler.setLevel( logging.FATAL)
            self.streamHandler.setLevel( logging.FATAL)
        else:
            self.fileHandler.setLevel( logging.INFO)
            self.streamHandler.setLevel( logging.INFO)
            

    @staticmethod
    def get_instance():
        if Logger.__instance == None:
            Logger.__instance = Logger()
            
        return Logger.__instance

        
    def debug(self, message):
        self.logg.debug(message)

    def info(self, message):
        self.logg.info(message)

    def warning(self, message):
        self.logg.warning(message)

    def error(self, message):
        self.logg.error(message + "\n---------------\n", exc_info=True)

    def critical(self, message):
        self.logg.critical(message + "\n---------------\n", exc_info=True)


