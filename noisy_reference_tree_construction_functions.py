import mysql.connector as mariadb
from ete3 import Tree, TreeNode
import pandas as pd
import pandas as pd
import polars as pl
import numpy as np
import fileinput
import random
import time
import os
import re

def nodes_indexation(PATH):
##This function takes as an input PATH to the nodes.dmp file
##and empty out list, and returns dataframe with parsed taxonomy levels
    df_nodes = []
    with open(PATH) as myfile:
        for line in myfile:
            line = line.split("\t")
            df_nodes.append([line[0], line[2], line[4]])
    #columns of taxonomy_df - tax_id, parent_id, tax_level
    df_nodes = pd.DataFrame(df_nodes, columns = ["tax_id", "parent_id", "tax_level"])
    df_nodes.set_index(["tax_id"], inplace = True)
    return df_nodes

def taxonomy_collector(df_taxonomy, initial_taxid, desired_level):
    ##This function takes as an input taxonomy levels dataframe,
    ##tax_id of initial group and desired level of subgroups
    ##and returns dataframe with tax_id, parent_id and tax_level of subgroups
    ## from given initial_taxid
    ## check that tax_level is valid
    available_taxonomy = pd.unique(df_taxonomy["tax_level"])
    print(available_taxonomy)
    if desired_level not in available_taxonomy:
        raise ValueError("desired_level must be one of %r." % available_taxonomy)
    initial_level = df_taxonomy.loc[initial_taxid]["tax_level"]
    desired_df = df_taxonomy[df_taxonomy["tax_level"] == desired_level]
    collected_taxonomy = []
    
    ##obtain all groups from desired level, and then go up 
    ##until root or initial_taxid achieved,
    ##return only groups fron initial_taxid
    for ind in desired_df.index:
        acc_initial = df_taxonomy.loc[str(ind)]
        acc = acc_initial
        while acc["tax_level"] != initial_level:
            acc = df_taxonomy.loc[str(acc["parent_id"])]
            if acc["parent_id"] == "1":
                break
        if acc.name == initial_taxid:
            collected_taxonomy.append([acc_initial.name, acc_initial["parent_id"], acc_initial["tax_level"]])
    #columns of collected_taxonomy_df - tax_id, parent_id, tax_level
    collected_taxonomy = pd.DataFrame(collected_taxonomy, columns = ["tax_id", "parent_id", "tax_level"])
    collected_taxonomy.set_index(["tax_id"], inplace = True)
    return collected_taxonomy

def tree_dataframe_construction(df_taxonomy, list_of_tax_id):
    #This function takes as an input df with taxonomy
    #and list with tax id and return dataframe
    #with columns based on upstream parent nodes

    #filter only necessary tax_ids
    desired_df = df_taxonomy.query("tax_id in @list_of_tax_id")
    index_counter = 0
    desired_df = desired_df.rename(columns={"parent_id": index_counter})
    #work up until all parens will be a root node
    while ~desired_df[index_counter].eq("1").all():
        new_df = df_taxonomy[df_taxonomy.index.isin(desired_df[index_counter])]
        map_df = new_df["parent_id"].to_dict()
        index_counter += 1
        desired_df[index_counter] = desired_df[index_counter-1].map(map_df)
    desired_df.drop(["tax_level"], axis = 1, inplace = True)
    print(desired_df)
    return desired_df

def recursive_tree(given_group, i, out_list):
    #This function is a recursive function
    #for tree construction based on dataframe
    #from tree_dataframe_construction
    if i != 1:
        #some recursion tricks, iterate for every group in nested groupby
        for group_name, group_data in given_group.groupby(i, as_index=False):
            out_group = recursive_tree(group_data, i - 1, [])
            out_list.append(out_group)
        return out_list
    else:
        #if group of level zero, return contents as list
        return given_group.index.tolist()
    
def flatten_tree(out_list):
    if isinstance(out_list, list):
        if len(out_list) == 1:
            return flatten_tree(out_list[0])
        else:
            return [flatten_tree(item) for item in out_list]
    else:
        return out_list
    
def get_table_names(mariadb_connection):
    #this function retreive all table names from
    #given swisspfam database
    cursor=mariadb_connection.cursor()
    cursor.execute("SHOW TABLES;")
    return [x[0] for x in cursor]

def taxid_to_mnemonics_species(PATH, nheader):
    ####This function reads speclist.txt file
    ##It takes as an input path to file with tax_id <-> mnemonic mapping
    ##Number of header lines to exclude from mapping file
    ##And dataframe to store the output data
    taxid_to_mnemonics_df = []
    with open(PATH, "r") as inp:
        for i in range(nheader):
            next(inp)
        for line in inp:
            line_n = line
            line = line.split(" ")
            #remove bad lines
            if len(line) != 1 and len([x for x in line if x.endswith(":")]) == 1:
                taxid_to_mnemonics_df.append([line[0], [x for x in line if x.endswith(":")][0][0:-1]])
        taxid_to_mnemonics_df = pd.DataFrame(taxid_to_mnemonics_df, columns = ["mnemonic", "tax_id"])
        taxid_to_mnemonics_df.set_index(["tax_id"], inplace = True)
    return taxid_to_mnemonics_df

def find_tax_level_for_tax_id(initial_taxid, input_df_taxonomy, desired_level):
    ##This function takes as an input df with taxonomical relations,
    ##taxid of entity and desired level of output
    ##and outputs the tax_id of the desired level
    if str(initial_taxid) in input_df_taxonomy.index:
        acc = input_df_taxonomy.loc[str(initial_taxid)]
        while acc["tax_level"] != desired_level:
            acc = input_df_taxonomy.loc[str(acc["parent_id"])]
            if acc["parent_id"] == "1":
                break
        if acc["tax_level"] == desired_level:
            return acc.name
    else:
        #print(f"{initial_taxid} not in df_taxonomy, return None")
        return None
    
def replace(match):
    value = int(match.group(0))
    print(str(map_df.get(str(value))))
    return str(map_df.get(str(value)))


def ete3_tree_from_df(df_with_relationships):
    ##This function constructs newick tree from given dataframe 
    ##with nodes relationships parsed from nodes.dmp 
    ##this function use tree structures from ete3 library
    
    nodes = {}
    root = TreeNode()
    root.name = '1'
    nodes['1'] = root

    for i, col in enumerate(reversed(df_with_relationships.columns)):
        for j, row in df_with_relationships.iterrows():
            node_label = str(row[col])
            if node_label in nodes:
                node = nodes[node_label]
            else:
                node = Tree()
                node.name = node_label
                nodes[node_label] = node
            if node.name != "1":
                child_node = node
                child_node.name = node_label
                nodes[node_label] = child_node
                parent_label = str(row[col+1])
                parent_node = nodes[parent_label]
                if child_node not in parent_node.children:
                    parent_node.add_child(child_node)
            
    newick_string = root.write(format=9)
    
    return newick_string