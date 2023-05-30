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

def filter_rf_output(PATH_TO_OUTPUT_DIR, remove = False):
    #This function read RF files in dir and filter them out
    #If RF value not calculated or fasta file unreachable
    if remove:
        print("Remove flag selected, following files will be deleted from the output:")
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
                    

def get_alignment_width_length(LIST_WITH_FASTA_PATHS):
    #This fucntion reads as an input a list of files
    #With paths for .fasta files (before or after denoising)
    #And returns Dataframe with width and length
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

def bootstrap_df_create(PATH_TO_BOOTSTRAP_RESULTS, dataframe):
    #This fucntion reads as an input the path to bootstrap .nwk
    #trees and calculate bootstrap statistics — amount of branches and 
    #sum of bootstrap values for a tree
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

def create_difference_dataframe(rf_df_denoised, rf_df_original):
    #This fucntion reads as an input 2 dataframes with denoised and original statistics
    #And returns one joint datarfame with total statistics
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
    rf_out_df["cutoff_percent"] = (rf_out_df["RF_original_length"] - rf_out_df["RF_denoised_length"])/rf_out_df["RF_original_length"]
    return rf_out_df


