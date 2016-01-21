# -*- coding: utf-8 -*-


##
# this class contains static methods aiming to simplify the usage of SQLalchemy results of query
#
class SQLAlchemyUtil(object):

    # format_list_column
    # ------------------
    #
    # Argument:
    #    - sql_alchemy_list : sql query result list (proceed on column).
    #
    # Loop on sql_alchemy_list and format the value, then append it into
    # the result list.
    #
    # Return result.
    #
    # @return list<String> or list<Integer> or list<Float>
    @staticmethod
    def format_list_column(sql_alchemy_list):
        result = []
        for e in sql_alchemy_list:
            result.append(e[0])

        return result

    # format_list_object
    # ------------------
    #
    # Argument:
    #    - sql_alchemy_list : sql query result list
    #                                            (proceed on complete object).
    #
    # Loop on sql_alchemy_list and format the value, then append it into
    # the result list.
    #
    # Return result.
    @staticmethod
    def format_list_object(sql_alchemy_list):
        result = []
        for e in sql_alchemy_list:
            result.append(e)

        return result
