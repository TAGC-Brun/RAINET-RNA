from _pyio import __metaclass__
from abc import ABCMeta, abstractmethod



class ExecutionStrategy( object):
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def execute(self):
        pass
    
    