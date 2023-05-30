import os
import re
import glob
import numpy as np
import polars as pl
import pandas as pd
import seaborn as sns
from scipy import stats
from Bio import AlignIO
import plotly.graph_objects as go
from scipy.stats import gaussian_kde
from plotly.subplots import make_subplots
from sklearn.neighbors import KernelDensity
from noisy_ds_construction_functions import *

#READ BOOTSTRAP LOGS
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

if __name__ == "__main__":

#CREATE LIST OF DATAFRAMES WITH RF
    for group in ["fungi", "metazoa", "viridiplantae", "proteobacteria", "actinobacteria", "archaea"]:
        print(group)
        PATH_TO_ORIGINAL_DIR = f"/data/ruslan_gumerov/{group}/{group}_results_original_bootstrap/"
        PATH_TO_DENOISED_DIR = f"/data/ruslan_gumerov/{group}/{group}_results_denoised_bootstrap/"
        PATH_TO_INITIAL_ALIGNMNENTS = f"/data/ruslan_gumerov/{group}/{group}_fasta_aligned/"

        #filter out bad files 
        if group == "archaea":
            filter_rf_output(PATH_TO_DENOISED_DIR, remove = True)
            filter_rf_output(PATH_TO_ORIGINAL_DIR, remove = True)
        else:
            filter_rf_output(PATH_TO_DENOISED_DIR, remove = False)
            filter_rf_output(PATH_TO_ORIGINAL_DIR, remove = False)

        #get RF vaLues
        rf_df_denoised = get_rf_from_log([x for x in  glob.glob(PATH_TO_DENOISED_DIR + "*tmp")])
        rf_df_original = get_rf_from_log([x for x in  glob.glob(PATH_TO_ORIGINAL_DIR + "*tmp")])

        #get alignment width and length
        fasta_df_denoised = get_alignment_width_length([x for x in glob.glob(PATH_TO_DENOISED_DIR + "*fas")])
        fasta_df_original = get_alignment_width_length([x for x in glob.glob(PATH_TO_INITIAL_ALIGNMNENTS + "*fasta")])

        #concat RF values and width/length of files
        rf_df_denoised = pd.concat([fasta_df_denoised, rf_df_denoised], axis = 1).dropna()
        rf_df_original = pd.concat([fasta_df_original, rf_df_original], axis = 1).dropna()

        #harmonize filenames
        rf_df_original = rf_df_original[rf_df_original.index.isin(rf_df_denoised.index)]
        rf_df_denoised = rf_df_denoised[rf_df_denoised.index.isin(rf_df_original.index)]

        
        #get domain info
        domains_table = pd.read_csv(f"/data/ruslan_gumerov/{group}/{group}_proteins_species_orders_mapping.csv", index_col=0)
        domains_table["apply"] = domains_table["apply"].map(clean_str)
        domains_table["taxa"] = domains_table[domains_table.columns[-1]].map(clean_str)
        PATH_TO_PROTEIN_DOMAIN_FILE = f"/data/ruslan_gumerov/{group}/{group}_domains_and_proteins_csv/"
        print(domains_table.apply(lambda row: get_domains_number(PATH_TO_PROTEIN_DOMAIN_FILE, row['apply'], row['taxa']),axis=1))
