import mysql.connector as mariadb
from Bio import AlignIO
import pandas as pd
import polars as pl
import numpy as np
import fileinput
import random
import time
import os

##MAIN
def nodes_indexation(PATH, df_nodes = []):
##This function takes as an input PATH to the nodes.dmp file
##and empty out list, and returns dataframe with parsed taxonomy levels
    with open(PATH) as myfile:
        for line in myfile:
            line = line.split("\t")
            df_nodes.append([line[0], line[2], line[4]])
    #columns of taxonomy_df - tax_id, parent_id, tax_level
    df_nodes = pd.DataFrame(df_nodes, columns = ["tax_id", "parent_id", "tax_level"])
    df_nodes.set_index(["tax_id"], inplace = True)
    return df_nodes

def names_indexation(PATH, df_names = []):
##This function takes as an input PATH to the names.dmp file
##and empty out list, and returns dataframe with tax_id and names
    with open(PATH) as myfile:
        for line in myfile:
            line = line.split("\t")
            df_names.append([line[0], line[2]])
    #columns of taxonomy_df - tax_id, name
    df_names = pd.DataFrame(df_names, columns = ["tax_id", "name"])
    return df_names

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

def find_tax_level_for_tax_id(initial_taxid, df_taxonomy, desired_level):
    ##This function takes as an input df with taxonomical relations,
    ##taxid of entity and desired level of output
    ##and outputs the tax_id of the desired level
    if str(initial_taxid) in df_taxonomy.index:
        acc = df_taxonomy.loc[str(initial_taxid)]
        while acc["tax_level"] != desired_level:
            acc = df_taxonomy.loc[str(acc["parent_id"])]
            if acc["parent_id"] == "1":
                break
        if acc["tax_level"] == desired_level:
            return acc.name
    else:
        print(f"{initial_taxid} not in df_taxonomy, return None")
        return None

def get_domains_list(mnemonic_list, connection, PATH_TO_DOMAIN_CSV):
    #This function takes as an input organism mnemonic and DB connection
    #Performs selection of ordered by domain_start groups of domains for each protein
    #And returns a list with ordered domains for this mnemonic
    cursor = connection.cursor()
    order_name = mnemonic_list.name
    mnemonic_list = mnemonic_list.values
    print(order_name)
    print(mnemonic_list)
    #if there is only one mnemonic in order
    if len(mnemonic_list) == 1:
        table_name = mnemonic_list[0]
        #get protein_code and domain sequence for mnemonic
        query = (f"SELECT protein_code, GROUP_CONCAT(domain_code ORDER BY domain_start) as domain_codes FROM {table_name} GROUP BY protein_code;")
        cursor.execute(query)
        #add mnemonic to protein_code
        q1 = [[element[0], element[1]] for element in cursor]
        proteins = [el[0] for el in q1]
        domains = [el[1] for el in q1]
        #we need only those domain architectures with exactly one protein
        mnemonic_out = []
        for protein, domain in zip(proteins, domains):
            if domains.count(domain) == 1:
                #add organism mnemonic
                mnemonic_out.append([[table_name + "_" + protein], [domain]])
        q1 = mnemonic_out
    else:
        q1 = []
        for mnemonic in mnemonic_list:
            table_name = mnemonic
            query = (f"SELECT protein_code, GROUP_CONCAT(domain_code ORDER BY domain_start) as domain_codes FROM {table_name} GROUP BY protein_code;")
            cursor.execute(query)
            query_out = [[element[0], element[1]] for element in cursor]
            proteins = [el[0] for el in query_out]
            domains = [el[1] for el in query_out]
            #we need only those domain architectures with exactly one protein
            query_out = []
            for protein, domain in zip(proteins, domains):
                if domains.count(domain) == 1:
                    #add organism mnemonic
                    query_out.append([[table_name + "_" + protein], [domain]])
            if not q1:
                #if common domain object for order is empty
                q1 = query_out
            else:
                #if not empty, add common protein codes for same domains
                domains_new = [element[1][0] for element in query_out]
                domains_old = [element[1][0] for element in q1]
                query_common = []
                for domain in set(domains_new).intersection(set(domains_old)):
                        modified_pair = [element for element in q1 if element[1][0] == domain][0]
                        new_protein_id = [element[0][0] for element in query_out if element[1][0] == domain]
                        modified_pair[0].append(new_protein_id[0])
                        query_common.append(modified_pair)
                for domain in set(domains_new).difference(set(domains_old)):
                        new_pair = [element for element in query_out if element[1][0] == domain][0]
                        query_common.append(new_pair)
                for domain in set(domains_old).difference(set(domains_new)):
                        old_pair = [element for element in q1 if element[1][0] == domain][0]
                        query_common.append(old_pair)
                q1 = query_common
    out_df = pd.DataFrame(q1, columns = ["proteins", "domains"])
    out_df.to_csv(f"{PATH_TO_DOMAIN_CSV}{order_name}.csv")
    cursor.close()
    return 0

def get_table_names(mariadb_connection):
    #this function retreive all table names from
    #given swisspfam database
    cursor=mariadb_connection.cursor()
    cursor.execute("SHOW TABLES;")
    return [x[0] for x in cursor]

def clean_str(s):
    #function for proper read of the csv columns from input
    return tuple(s.strip("[]").replace("'","").replace(" ","").split(","))

def selector(x):
    #function for random selection of proteins from input
    return pl.Series([random.choice(proteins) for proteins in x[0]])

def filter_rf_output(PATH_TO_OUTPUT_DIR, remove = False):
    if remove:
        print("Remove flag selected, following files will be deleted from the output:")
    #This function read file in dir and filter them out
    #If RF value not calculated or fasta file unreachable
    for tmp_file in os.listdir(PATH_TO_OUTPUT_DIR):
        if os.path.splitext(tmp_file)[1] == ".tmp":
            with open(f"{PATH_TO_OUTPUT_DIR}{tmp_file}", 'r') as f:
                #Try to read RF 
                try:
                    f.readline()
                    rf_dist = float(f.readline().strip())
                except Exception:
                    #If RF value not in the first line, it may also be in the second one
                    with open(f"{PATH_TO_OUTPUT_DIR}{tmp_file}", 'r') as f:
                        try:
                            rf_dist = float(f.readline().strip())
                        except Exception:
                            #If RF value not in the second line, exclude file from the analysis
                            print(f"problems with file {tmp_file}")
                            if remove:
                                os.remove(f"{PATH_TO_OUTPUT_DIR}{tmp_file}")
    #Try to read fasta file
    for fasta_file in os.listdir(PATH_TO_OUTPUT_DIR):
        if os.path.splitext(fasta_file)[1] == ".fas":
            try:
                alignment = AlignIO.read(f"{PATH_TO_OUTPUT_DIR}{fasta_file}", "fasta")
            #If not readable remove from the analysis
            except ValueError:
                print(f"problems with file {fasta_file}")
                if remove:
                    os.remove(f"{PATH_TO_OUTPUT_DIR}{fasta_file}")

def get_rf_from_log(LIST_WITH_RF_PATHS):
    rf_dict = dict()
    rf_dict["name"] = ['_'.join(os.path.basename(x).split("_")[0:2]) for x in LIST_WITH_RF_PATHS]
    rf_list = []
    for tmp_file in LIST_WITH_RF_PATHS:
        with open(tmp_file, 'r') as f:
            try:
                f.readline()
                rf_list.append(float(f.readline().strip()))
            except Exception:
                with open(tmp_file, 'r') as f:
                    rf_list.append(float(f.readline().strip()))
    rf_dict["RF"] = rf_list
    return pd.DataFrame.from_dict(rf_dict).set_index("name")

def get_alignment_width_length(LIST_WITH_FASTA_PATHS):
    fasta_dict = dict()
    fasta_dict["name"] = ['_'.join(os.path.basename(x).replace(".", "_").split("_")[0:2]) for x in LIST_WITH_FASTA_PATHS]
    len_list = []
    wid_list = []
    for fasta_file in LIST_WITH_FASTA_PATHS:
        alignment = AlignIO.read(fasta_file, "fasta")
        len_list.append(alignment.get_alignment_length())
        wid_list.append(len(alignment))
    fasta_dict["sequence_length"] = len_list
    fasta_dict["sequence_width"] = wid_list
    return pd.DataFrame.from_dict(fasta_dict).set_index("name")


def get_domains_number(PATH_TO_PROTEIN_DOMAIN_FILE, protein_name_list, group_tax_id_list):
    #apply this function to domains_and_proteins.csv file to get exact amount of domains
    return_list = []
    for protein_name, group_tax_id in zip(protein_name_list, group_tax_id_list):
        file = pd.read_csv(f"{PATH_TO_PROTEIN_DOMAIN_FILE}{group_tax_id}.csv", index_col=0)
        file["proteins"] = file["proteins"]
        file["domains"] = file["domains"]
        print(file["proteins"])
        num_domains = file[file["proteins"].str.contains(protein_name)]["domains"].str.split(",").str.len().values[0]
        return_list.append(num_domains)
    
    return set(return_list)

def create_difference_dataframe(rf_df_denoised, rf_df_original):
    rf_difference_column = rf_df_original.sort_index()["RF"] - rf_df_denoised.sort_index()["RF"]
    print(f"0 difference: {sum(rf_difference_column == 0)}, "
          f"original RF higher than denoised: {sum(rf_difference_column > 0)}, ", 
          f"original RF lower than denoised: {sum(rf_difference_column < 0)}")
    rf_out_df = pd.DataFrame({"RF_original_length" : rf_df_original["sequence_length"], 
                              "RF_denoised_length" : rf_df_denoised["sequence_length"], 
                              "RF_original_width" : rf_df_original["sequence_width"], 
                              "RF_original" : rf_df_original["RF"], 
                              "RF_denoised" : rf_df_denoised["RF"], 
                              "RF_difference" : rf_difference_column,
                              "RF_bool" : None})
    rf_out_df.loc[rf_out_df["RF_difference"] == 0, "RF_bool"] = "zero"
    rf_out_df.loc[rf_out_df["RF_difference"] > 0, "RF_bool"] = "original_higher"
    rf_out_df.loc[rf_out_df["RF_difference"] < 0, "RF_bool"] = "original_lower"
    rf_out_df["cutoff_percent"] = rf_out_df["RF_original_length"] / rf_out_df["RF_denoised_length"]
    return rf_out_df

def bootstrap_df_create(PATH_TO_BOOTSTRAP_RESULTS, dataframe):
    file_list = []
    branch_list = []
    bootstrap_list = []
    for file in os.listdir(PATH_TO_BOOTSTRAP_RESULTS):
        if file.endswith("tree.nwk"):
            with open(f"{PATH_TO_BOOTSTRAP_RESULTS}{file}", "r") as inp:
                first_line = inp.readline()
                matches = re.findall(r'\)(\d+):', first_line)
                file_list.append(file)
                branch_list.append(len(matches))
                bootstrap_list.append(sum(list(map(int, matches))))
                
    dataframe["file"] = file_list
    dataframe["branch_length"] = branch_list
    dataframe["bootstrap"] = bootstrap_list
    return dataframe
