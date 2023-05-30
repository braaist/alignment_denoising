#get main df
from ../noisy_ds_construction_functions import *
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
OUT_FOLDER_PATH = "/data/ruslan_gumerov/metazoa/"
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
    print("STEP 1. File indexation and taxonomy_df construction.")
    df_nodes_with_parents = nodes_indexation(PATH_nodes)
    #get dataframe tax_id <-> mnemonics
    taxid_to_mnemonics_df = taxid_to_mnemonics_species(PATH_speclist, 59)
    #get dataframe tax_id <-> parent_tax_id <-> tax_level
    #for all species in parent_tax_level
    mariadb_connection = mariadb.connect(user="ruslan_gumerov",
                                     password="",
                                     database="swisspfam",
                                     unix_socket="/run/mysqld/mysqld.sock")
    swisspfam_mnemonics = get_table_names(mariadb_connection)
    taxid_to_mnemonics_df.query("mnemonic in @swisspfam_mnemonics", inplace = True)
    taxid_to_mnemonics_df[parental_tax_level] = taxid_to_mnemonics_df.index.to_series().apply(find_tax_level_for_tax_id,
                                                                          df_taxonomy = df_nodes_with_parents, 
                                                                          desired_level = parental_tax_level)
    taxid_to_mnemonics_df = taxid_to_mnemonics_df[taxid_to_mnemonics_df[parental_tax_level] == parental_tax_id]
    taxid_to_mnemonics_df[aggregation_tax_level] = taxid_to_mnemonics_df.index.to_series().apply(find_tax_level_for_tax_id,
                                                                          df_taxonomy = df_nodes_with_parents, 
                                                                          desired_level = aggregation_tax_level)
    taxid_to_mnemonics_df.dropna(inplace=True)

    taxid_to_mnemonics_df.groupby(aggregation_tax_level)["mnemonic"].apply(get_domains_list,  mariadb_connection, PATH_TO_DOMAIN_CSV)
    mariadb_connection.close()


    print("STEP 2. Filter domains and get random proteins for a given group.")
    random.seed(5051)
    #read domains files
    domains_df = pl.read_csv(f"{PATH_TO_DOMAIN_CSV}*.csv")

    #read df with results
    result = domains_df.select([
        pl.col("proteins").apply(clean_str),
        pl.col("domains").apply(clean_str)
    ]).with_columns(pl.lit(1).alias('amount')).groupby("domains").agg([pl.col("proteins"), pl.col("amount").sum()])

    #filter only domains with amount >= 15
    result_filtered = result.filter(pl.col("amount") >= 15)

    #get pandas series of lists with selected proteins
    proteins = result_filtered.select([
        pl.col("proteins")
    ])
    proteins_out = proteins.apply(selector).to_pandas()

    #read dataframe with mnemonics <-> orders
    #and create mapping dictionary
    map_df = taxid_to_mnemonics_df.drop([parental_tax_level], axis = 1).set_index('mnemonic')[aggregation_tax_level].to_dict()

    #add columns with corresponding species
    # and orders list
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
