# Class Table gets as input a table and allows to perform several functions for parsing it.

import operator
from core.util.log.Logger import Logger
from core.util.exception.ParsingFileException import ParsingFileException



class TableWrapper(): 
    
    # get_column
    # --------------------------------------------------
    #
    # this method extracts one particular column (index) 
    #
    # Arguments:
    #    - table: 
    #    - index: column index that you want get
    #    - start and end: indexes that can be specify in order to get onw only part of column
    #
    # Returns:
    #    - column: list
    @staticmethod
    def get_column(table, index, start=None, end=None):
        
        # Indexes check 
        if index > len(table[0]):
            Logger.get_instance().error(" TableWrapper.get_column : the column index is greater then column number of table\n")
            raise ParsingFileException(" TableWrapper.get_column : the column index is greater then column number of table\n")
        if start != None and end != None:
            if start == end:
                Logger.get_instance().error(" TableWrapper.get_column : start and end indexes can't be equal\n")
                raise ParsingFileException(" TableWrapper.get_column :  start and end indexes can't be equal\n")
            elif start > end:
                Logger.get_instance().error(" TableWrapper.get_column : start index can't be greater than end index\n")
                raise ParsingFileException(" TableWrapper.get_column : start index can't be greater than end index\n")
            else:
                Logger.get_instance().info(" TableWrapper.get_column : start and end indexes are correct\n" +'start: '+str(start)+ ', end: '+str(end))
        elif start != None and end == None:
            Logger.get_instance().error(" TableWrapper.get_column : start index is greater than table length\n")
            raise ParsingFileException(" TableWrapper.get_column : start index is greater than table length\n")

        # in according to combination of index, start and end indexes the method returns 
        # a specific column
        if start == None and end == None:
            columns = zip(*table)
            return list(columns[index])
        elif start != None and end != None:
            columns = zip(*table)
            return list(columns[index][start:end])
        elif start != None and end == None:
            columns = zip(*table)
            return list(columns[index][start:])
        elif start == None and end != None:
            columns = zip(*table)
            return list(columns[index][0:end])

    # del_column
    # --------------------------------------------
    #
    # This method deletes a column of table
    #
    # Arguments:
    #    - table
    #    - index: column to be delete
    #
    # Returns: 
    #    - New table without one column
    @staticmethod
    def del_column(table, index):
        if index > len(table[0]):
            Logger.get_instance().error(" TableWrapper.get_column : the column index is greater then column number of table\n")
            raise ParsingFileException(" TableWrapper.get_column : the column index is greater then column number of table\n")
        new_table = []
        for row in table:
            row.pop(index)
            new_table.append(row)
        return new_table
    
    
    # inv_column
    # ----------------------------------------------
    #
    # the inv_column get input two indices and inverts the respective columns
    #
    # Arguments:
    #    - table
    #    - index1 and index2: the indexes of column that you want 
    # Returns:
    #    - new_table 
    @staticmethod
    def inv_column(table, index1, index2):
        if index1 > len(table[0]) or index2 > len(table[0]):
            Logger.get_instance().error(" TableWrapper.get_column : the column index is greater then column number of table\n")
            raise ParsingFileException(" TableWrapper.get_column : the column index is greater then column number of table\n")
        new_columns = []
        new_table = []
        num = len(table[0])
        for i in range(num):
            if i != index1 and i != index2:                
                new_columns.append(TableWrapper.get_column(table,i))
            elif i == index1:
                new_columns.append(TableWrapper.get_column(table, index2))
            elif i == index2:
                new_columns.append(TableWrapper.get_column(table, index1))
        for item in zip(*new_columns):
            item = list(item)
            new_table.append(item)
        return new_table
    
    # add_column
    # ------------------------------------------------------
    #
    # this method add one columns before first column
    #
    # Arguments:
    #    - table: 
    #    - column: list of item that you want add at the table
    # Returns:
    #    - New table
    @staticmethod
    def add_column(table, column):
        new_table = []
        for i, item in enumerate(column):
            row=[]
            row.append(item)
            for f in range(len(table[i])):                
                row.append(table[i][f])
            new_table.append(row)
        return new_table 
    
    #
    # add_end_column
    # ------------------------------------------------------
    #
    # this method add one column after last column 
    #
    # Arguments:
    #    - table: 
    #    - column: list of item that you want add at the table
    # Returns:
    #    - New table
    @staticmethod
    def add_end_column(table, column):
        new_table = []
        for i, item in enumerate(column):
            row=[]
            for f in range(len(table[i])):                
                row.append(table[i][f])
            # the items of new column are added at the end of each row
            row.append(item)
            new_table.append(row)
        return new_table 

    #
    # make_dictionary
    # --------------------------------------------------------
    #
    # This method makes a dictionary 
    # the keys are the terms of first column of table and the value are those of second column 
    # If the first column present the same item the value are appended to list values of key
    #
    # Arguments: 
    #     - table: table with almost 2 columns
    #     - term: number of values to be add
    # Returns:
    #     - dictionary = {keys: [item1, item2..itemN]}         
    @staticmethod
    def make_dictionary(table, term=1):
        dictionary={}
        for row in table:

            # term = 1
            # Table Example: [[keyA, value1],[keyB,value1],[keyC,value1],[keyA,value2]]
            # Note in this case can be duplicate keys
            
            if term==1:         
                if row[0] in dictionary:
                    dictionary[row[0]].append(row[1])
                else:
                    dictionary[row[0]]=[row[1]]
            # term > 1
            # Table Example: [[keyA, value1, feature1],[keyB,value1,feature1],[keyC,value1,feature1]]
            # Note: in this case there aren't duplicate keys
            else:
                for item in row[1:term+1]:
                    if row[0] in dictionary:
                        dictionary[row[0]].append(item)
                    else:
                        dictionary[row[0]]=[item]
        return dictionary

    #this method returns the table (provided as input) sorted by column pass to argument "col"
    #The default column of sorting is the first one (index = 0)
    #
    @staticmethod
    def sort_table(table, col=0): 
        return sorted(table, key=operator.itemgetter(col), reverse=True)

