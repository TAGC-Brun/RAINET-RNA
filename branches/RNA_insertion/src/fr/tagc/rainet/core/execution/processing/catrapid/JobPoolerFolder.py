
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
    
    #bash runlibrary.both.sh prot/prot.lib rna/rna.lib 0 >log.out 2>log.err
    
    cmd = "cd %s; %s" % ( argList[0], argList[1]) 
    # TODO: use subprocess? get logging , perhaps use my function
    os.system( cmd)
    return cmd
 
class JobPoolerFolder( object):
    
    def __init__(self, input_folder, command_string, num_threads, output_folder):

        self.inputFolder = input_folder
        
        self.commandString = command_string

        # self.scriptPath = script_path

        self.numThreads = num_threads
        
        self.outputFolder = output_folder

    # #
    # Function spawning jobs in a pool
    def parallel( self):
                
        # get list of folders to be processed
        foldersToProcess = sorted( glob.glob( self.inputFolder + "/*") )
        
        # add script path before each folder, separeted by "|"
        runList = [ "||".join( [folder, self.commandString] ) for folder in foldersToProcess]
        
        # initialise pool    
        pool = Pool( self.numThreads)
            
        # note that second argument needs to be iterable, and a job will be launched for each        
        results = pool.map( run_command_list_args, runList)
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
    parser.add_argument('outputFolder', metavar='outputFolder', type=str, help='Folder where to write output files.')
    
    #gets the arguments
    args = parser.parse_args( ) 

    pooler = JobPoolerFolder( args.inputFolder, args.commandString, args.numThreads, args.outputFolder)

    results = pooler.parallel( )
    
#     for result in results:
#         print result


    print "FINISHED!"
