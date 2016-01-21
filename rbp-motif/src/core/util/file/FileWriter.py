# This class writes into file a column or a table

from core.util.file.FileUtils import FileUtils
from core.util.exception.RbpmotifException import RbpmotifException

class FileWriter(): 
    
    # write_table
    # ------------------------------------------
    # write a list object in a file
    #
    # Arguments:
    #    - file: file path 
    #    - obj : list object (simple list or table)
    #    - symbol : symbol that separates the table items 
    # 
    # Return: 
    #    - file 
    @staticmethod
    def write_table(file_path, obj, symbol='\t'):
        outfile = FileUtils.open_text_w(file_path)
        # It checks if the item of list are string, integer or float: if the items are ones of these type 
        # the items are written in a file in column format
        # just in the case the item of list are again a list (that is a table object) the object will be
        # written in a table format where the columns are separate by a symbol that can be chosen as input argument
        try: 
            for item in obj:
                if type(item) != str and type(item) != int and type(item) != float:
                    length = len(item)
                    for i in range(length):
                        outfile.write(str(item[i]) + symbol)
                    outfile.write('\n')
                elif type(item) == str:
                    outfile.write(str(item) + '\n')
                elif type(item) == int or type(item) == float:
                    outfile.write(str(item) + '\n')
        except:
            raise RbpmotifException("FileWriter.write_table: Error in file writing :" + file_path)