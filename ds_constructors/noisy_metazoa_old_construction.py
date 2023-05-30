#get main df
from noisy_ds_construction_functions import *
import mysql.connector as mariadb
import pandas as pd
import polars as pl
import numpy as np
import fileinput
import random
import time
import os

#tax level and group name 
#of group for main level (kingdom)
GROUP_NAME = "METAZOA"
parental_tax_level = "kingdom"
parental_tax_id = "33208"
#initial tax level for dataset creation (species)
initial_tax_level = "species"
#tax level for proximity penalisation (genus)
aggregation_tax_level = "order"
#path to nodes.dmp file
PATH_nodes = "/data/DBs/Taxonomy/13-05-2022/nodes.dmp"
#path to speclist file
PATH_speclist = "/data/DBs/Uniprot/speclist.txt"
#path to main output folder
OUT_FOLDER_PATH = "/data/ruslan_gumerov/metazoa_old_test/"
#path to folder with domains and proteins for given group_name and tax_levels
#path to output fasta folder
PATH_TO_OUTPUT_FASTA_DIR = OUT_FOLDER_PATH + "metazoa_fasta/"
#path to the table with mapping between randomly selected proteins, their species and genuses
PATH_TO_MAIN_SELECTION_TABLE = OUT_FOLDER_PATH + "metazoa_proteins_species_orders_mapping.csv"
#path to the output .csv files with domains and proteins for given genus
PATH_TO_DOMAIN_CSV = OUT_FOLDER_PATH + "metazoa_domains_and_proteins_csv/"

#path creation
for path in [OUT_FOLDER_PATH, PATH_TO_OUTPUT_FASTA_DIR, PATH_TO_DOMAIN_CSV]:
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Path created: {path}")

if __name__ == "__main__":
    #get dataframe tax_id <-> parent_tax_id <-> parent_tax_level
    #print("STEP 1. File indexation and taxonomy_df construction.")
    #df_nodes_with_parents = nodes_indexation(PATH_nodes)
    #get dataframe tax_id <-> mnemonics
    #taxid_to_mnemonics_df = taxid_to_mnemonics_species(PATH_speclist, 59)
    #get dataframe tax_id <-> parent_tax_id <-> tax_level
    #for all species in parent_tax_level
    #mariadb_connection = mariadb.connect(user="ruslan_gumerov",
    #                                 password="",
    #                                 database="swisspfam",
    #                                 unix_socket="/run/mysqld/mysqld.sock")
    #swisspfam_mnemonics = get_table_names(mariadb_connection)
    #taxid_to_mnemonics_df.query("mnemonic in @swisspfam_mnemonics", inplace = True)
    #taxid_to_mnemonics_df[parental_tax_level] = taxid_to_mnemonics_df.index.to_series().apply(find_tax_level_for_tax_id,
    #                                                                      df_taxonomy = df_nodes_with_parents, 
    #                                                                      desired_level = parental_tax_level)
    #taxid_to_mnemonics_df = taxid_to_mnemonics_df[taxid_to_mnemonics_df[parental_tax_level] == parental_tax_id]
    #taxid_to_mnemonics_df[aggregation_tax_level] = taxid_to_mnemonics_df.index.to_series().apply(find_tax_level_for_tax_id,
                                                                          df_taxonomy = df_nodes_with_parents, 
                                                                          desired_level = aggregation_tax_level)
    #taxid_to_mnemonics_df.dropna(inplace=True)

    #taxid_to_mnemonics_df.groupby(aggregation_tax_level)["mnemonic"].apply(get_domains_list,  mariadb_connection, PATH_TO_DOMAIN_CSV)
    #mariadb_connection.close()


    print("STEP 2. Get old domains.")

    Metazoa_80 = pd.read_csv("/data/ruslan_gumerov/Metazoa_old.tsv", sep = "\t")
    Metazoa_80.columns = Metazoa_80.columns[1:].insert(0, "domain_group")
    Metazoa_80["domain_name"] = Metazoa_80["sp"] + "_" +Metazoa_80["prot"]
    Metazoa_80 = Metazoa_80.groupby(['domain_group'])['domain_name'].apply(','.join).reset_index()
    proteins_out = pd.DataFrame({"apply" : Metazoa_80["domain_name"].str.split(",")})

    #read mnemonics file
    taxid_to_mnemonics_df = pd.read_csv("/data/ruslan_gumerov/mnemonics_df_metazoa_test.csv")
    map_df = taxid_to_mnemonics_df.drop([parental_tax_level], axis = 1).set_index('mnemonic')[aggregation_tax_level].to_dict()
    
    #exclude different mnemonics
    out_set = set()
    for x in proteins_out["apply"]:
        for y in x:
            if y.split("_")[0] not in map_df:
                out_set.add(y.split("_")[0])
    item_list = []
    for item in proteins_out["apply"]:
        protein_list = []
        for protein in item:
            if not any(protein.startswith(p) for p in out_set):
                protein_list.append(protein)
        item_list.append(protein_list)

    proteins_out["apply"] = item_list
    proteins_out["species"] = [[y.split("_")[0] for y in x] for x in proteins_out["apply"]]
    proteins_out[aggregation_tax_level] = [[map_df[y.split("_")[0]] for y in x] for x in proteins_out["apply"]]


    proteins_out.to_csv(PATH_TO_MAIN_SELECTION_TABLE)

    print("STEP 3. Create .fasta files for given proteins.")
    #main .fasta creation loop
    #iterates through orders and proteins from
    #proteins_out df and run seqret
    counter_fasta = 0
    for proteins, orders in zip(proteins_out["apply"], proteins_out[aggregation_tax_level]):
        proteins = [x.split("_")[1] + "_" + x.split("_")[0] for x in proteins]
        mapping_species_orders = {">"+proteins[i]: orders[i] for i in range(len(proteins))}
        os.system(f"seqret pfamseq:{','.join(proteins)} -out {PATH_TO_OUTPUT_FASTA_DIR}/{counter_fasta}.fasta")
        #open created fasta and replace protein codes by orders
        for line in fileinput.input(f"{PATH_TO_OUTPUT_FASTA_DIR}/{counter_fasta}.fasta", inplace=True):
            if line.startswith(">"):
                print(">" + str(mapping_species_orders[line.split(" ")[0]]), end = '\n')
            else:
                print(line, end = '')
        counter_fasta += 1
