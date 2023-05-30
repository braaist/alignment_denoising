from NOISY_reference_tree_construction_functions import *
import mysql.connector as mariadb
import pandas as pd
PREFIX = "/data/ruslan_gumerov/"

#Creation of Pfam data table
if __name__ == "__main__":
    PATH_nodes = "/data/DBs/Taxonomy/13-05-2022/nodes.dmp"
    PATH_speclist = "/data/DBs/Uniprot/speclist.txt"
    
    df_initial = nodes_indexation(PATH_nodes)
    taxid_to_mnemonics_df = taxid_to_mnemonics_species(PATH_speclist, 59)
    mnemonics = list(set(taxid_to_mnemonics_df["mnemonic"]))
    
    #convert mnemonics to tax_ids
    swisspfam_id = []
    no_taxid = []
    for x in mnemonics:
        if not taxid_to_mnemonics_df[taxid_to_mnemonics_df["mnemonic"] == x].empty:
            swisspfam_id.append(taxid_to_mnemonics_df[taxid_to_mnemonics_df["mnemonic"] == x].index[0])
        else:
            swisspfam_id.append(None)
            
    # create and save final DF
    mnemonic_taxid_df = pd.DataFrame({"swisspfam_mnemonic" : mnemonics,
                           "tax_id" : swisspfam_id})

    mnemonic_taxid_df.to_csv(f"{PREFIX}noisy/mnemonic_taxid_df.csv")