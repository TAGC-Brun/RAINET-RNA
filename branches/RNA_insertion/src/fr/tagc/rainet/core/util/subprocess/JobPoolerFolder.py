
import os
import glob
import argparse

from multiprocessing import Pool

#===============================================================================
# Started 13-May-2016 
# Diogo Ribeiro
DESC_COMMENT = "Script to run jobs in a thread pool. Processes wanted command for each folder in input folder."
#===============================================================================

# #
# Function to be run for each parallel job. 
# This method must be unbound (e.g. not inside a class)
# Normally only one argument supported, here we split a string to harbour several arguments
# @param argString : string containing several arguments separated by "||"
def run_command_list_args( argString):
    
    argList = argString.split( "||")
        
    cmd = "sleep %s; cd %s; %s" % ( argList[0], argList[1], argList[2]) 
    # TODO: use subprocess? SubprocessUtil.run_command( command)

    os.system( cmd)

    return cmd
 
class JobPoolerFolder( object):
    
    def __init__(self, input_folder, command_string, num_threads, sleep_time):

        self.inputFolder = input_folder
        
        self.commandString = command_string

        self.numThreads = num_threads
        
        self.sleepTime = sleep_time
        
    # #
    # Function spawning jobs in a pool
    def parallel( self):
                
        # get list of folders to be processed
        foldersToProcess =  sorted( glob.glob( self.inputFolder + "/*"))
        
        # create sleep times list
        # purpose is to have the first jobs start at different times, so that they get initially defased.
        sleepTimes = {}
        for (i, folder) in enumerate( foldersToProcess ):
            if i < self.numThreads:
                sleepTimes[ folder] = str( i * self.sleepTime)
            else:
                sleepTimes[ folder] = str( 0)

        # create list of jobs, add arguments separated by "||"
        runList = [ "||".join( [ sleepTimes[ folder], folder, self.commandString] ) for folder in foldersToProcess]
        
        # initialise pool    
        pool = Pool( self.numThreads)
             
        # note that second argument needs to be iterable, and a job will be launched for each        
        results = pool.map( run_command_list_args, runList, chunksize = 1)
        pool.close()
        pool.join()
         
        return results
 
if __name__ == "__main__":
    
    print "STARTING!"
    
    #===============================================================================
    # Get input arguments, initialise class
    #===============================================================================
    parser = argparse.ArgumentParser(description= DESC_COMMENT) 

    # positional args
    parser.add_argument('inputFolder', metavar='inputFolder', type=str,
                         help='Folder containing the folders to be processed.')
    parser.add_argument('commandString', metavar='commandString', type=str,
                         help='Command to be run for each folder.')
#     parser.add_argument('scriptPath', metavar='scriptPath', type=str,
#                          help='Script to be executed for each folder in inputFolder. Needs to be executable.')
    parser.add_argument('numThreads', metavar='numThreads', type=int,
                         help='Number of threads to use.')
#     parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
    parser.add_argument('--sleepTime', metavar='sleepTime', type=int, default = 120,
                         help='Sleep timer for the start jobs, so that jobs get dephased from each other.')
    
    #gets the arguments
    args = parser.parse_args( ) 

    pooler = JobPoolerFolder( args.inputFolder, args.commandString, args.numThreads, args.sleepTime)

    results = pooler.parallel( )
    
#     for result in results:
#         print result


    print "FINISHED!"
