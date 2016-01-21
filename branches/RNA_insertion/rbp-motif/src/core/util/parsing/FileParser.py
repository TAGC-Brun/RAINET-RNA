# This class gets as input a file that you want to parse

from core.util.file.FileUtils import FileUtils

class FileParser():
    
    #
    # read_file
    # -----------------------------------------------
    #
    # this method reads a file and put all rows into list
    #
    # Argument:
    #    -file name
    # Returns:
    #    -listfile: list of file rows
    #
    @staticmethod
    def read_file(namefile):
        f = FileUtils.open_text_r(namefile)
        listfile = []
        for line in f:
            item = line.strip()
            listfile.append(item)
        return listfile
    
    #
    # make_table
    # ----------------------------------------------------------
    #
    # this method splits the table according to symbol provided as input
    #
    # Arguments:
    #    - namefile
    #    - symbol: 
    #    - skip:  allows to skip one or more header lines
    # Returns:
    #    - table object of file read
    @staticmethod
    def make_table(namefile, symbol='\t', skip=None):
        list_file = FileParser.read_file(namefile)
        table = []
        if skip==None:
            ind=0
        elif skip != None:
            ind = skip
        for item in list_file[ind:]:
            row = item.split(symbol)
            table.append(row)
        return table
    
    # 
    # merge_file
    # ----------------------------------------------------
    #
    # this method merges two files
    #
    @staticmethod
    def merge_file(namefile_1,namefile_2, new_namefile):
        file1 = FileUtils.open_text_r(namefile_1)
        file2 = FileUtils.open_text_r(namefile_2)
        new_file = FileUtils.open_text_a(new_namefile)
        text_1 = file1.read()
        text_2 = file2.read()
        new_file.write(text_1)
        new_file.write(text_2)
        new_file.close()







