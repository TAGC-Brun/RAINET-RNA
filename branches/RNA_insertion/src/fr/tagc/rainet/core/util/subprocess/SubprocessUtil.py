
import subprocess

from fr.tagc.rainet.core.util.exception.RainetException import RainetException
from fr.tagc.rainet.core.util.log.Logger import Logger

# #
# This class contains utilities about regular expression and pattern
class SubprocessUtil( object ):    

    # #
    # Method to run external commands and retrieve standard output and error.
    # Sends command and waits for it to finish. Then logs STDERR as a Logger warning,
    # STDOUT as Logger info and return code != 0 as Logger error.
    # 
    # @param command : string - the bash-like command to be ran
    @staticmethod  
    def run_command(self, command, decoder = "UTF-8"):

        try:
            command = str( command)
        except TypeError:
            raise RainetException( )

        Logger.get_instance().info( "\nrun_command : Starting command : "+command)
        
        p = subprocess.Popen( command, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        
        while True:
            if p.poll() != None:
                break

        stderrText = p.stderr.read().decode( decoder).strip()
        stdoutText = p.stdout.read().decode( decoder).strip()
    
        commentLine = "\n##################\n"
    
        if len(stderrText) > 0:
            Logger.get_instance().warning( "\nSubprocessUtil.run_command :\n%s# %s%s%s%s" % ( commentLine, "STDERR", commentLine, stderrText, commentLine) )

        if p.returncode != 0:
            Logger.get_instance().error( "\nSubprocessUtil.run_command :\n%s# %s%s%s%s" % ( commentLine, "ERROR (return code != 0)", commentLine, str(p.returncode), commentLine) )
        
        if len(stdoutText) > 0:
            Logger.get_instance().info( "\nSubprocessUtil.run_command :\n%s# %s%s%s%s" % ( commentLine, "STDOUT", commentLine, stdoutText, commentLine) )

        return None

