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

def noisy_density_before_after(rf_df_original, rf_df_denoised, name = "Not specified"):
    tmp = pd.DataFrame({ "log_length" : np.log(rf_df_original["sequence_length"])})
    sns_plot1 = sns.kdeplot(data=tmp, x="log_length")
    tmp = pd.DataFrame({ "log_length" : np.log(rf_df_denoised["sequence_length"])})
    sns_plot2 = sns.kdeplot(data=tmp, x="log_length")
    del tmp
    
    density_fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    density_fig.add_trace(
        go.Histogram(
            x=np.log(rf_df_denoised["sequence_length"]),
            name='после noisy', # name used in legend and hover labels
            marker_color='#EF553B',
            nbinsx=50,
            opacity=0.75),
        secondary_y=False)
    
    density_fig.add_trace(
        go.Histogram(
            x=np.log(rf_df_original["sequence_length"]),
            name='исходные',
            marker_color='#00CC96',
            nbinsx=50,
            opacity=0.75),
        secondary_y=False)
    
    density_fig.add_trace(
        go.Scatter(
            x=sns_plot1.get_lines()[0].get_data()[0], 
            y=sns_plot1.get_lines()[0].get_data()[1], 
            mode="lines", 
            marker_color='#00CC96',
            name="исходные",
            line = dict(width = 1.5, dash='dash')),
        secondary_y=True)
    
    density_fig.add_trace(
        go.Scatter(
            x=sns_plot2.get_lines()[1].get_data()[0], 
            y=sns_plot2.get_lines()[1].get_data()[1], 
            mode="lines", 
            marker_color='#EF553B',
            name="после noisy",
            line = dict(width = 1.5, dash='dash')),
        secondary_y=True)

    density_fig.update_layout(
        title= f"Плотность длин для разных классов в {name}",
        xaxis_title="Длина последовательности (log)",
        template="none",
    )
    
    density_fig.update_yaxes(title_text="Количество последовательностей", rangemode='tozero', secondary_y=False)
    density_fig.update_yaxes(title_text="Плотность", secondary_y=True, showgrid=False,
                     visible=False, showticklabels=False, rangemode='tozero', tickmode="auto")

    del sns_plot1
    del sns_plot2

    return density_fig

def noisy_cutoff_density(rf_df_original, rf_df_denoised, name = "Not specified"):
    tmp = pd.DataFrame({ "cutoff" : rf_df_original["sequence_length"] - rf_df_denoised["sequence_length"]})
    sns_plot = sns.kdeplot(data=tmp, x="cutoff")
    del tmp
    
    cutoff_fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    cutoff_fig.add_trace(
        go.Histogram(
            x=rf_df_original["sequence_length"] - rf_df_denoised["sequence_length"],
            name='original',
            marker_color='#636EFA',
            nbinsx=50,
            opacity=0.75),
        secondary_y=False)

    cutoff_fig.update_traces(marker_line_width=1,marker_line_color="black")

    cutoff_fig.update_layout(
        title= f"Распределение количества колонок, отрезанных noisy в {name}",
        yaxis_title="Количество последовательностей",
        xaxis_title="Размер отрезанной части",
        template="none",
    )
    
    cutoff_fig.update_yaxes(title_text="Cutoff size", rangemode='tozero', secondary_y=False)
    cutoff_fig.update_yaxes(title_text="Cutoff size", secondary_y=True, showgrid=False,
                     visible=False, showticklabels=False, rangemode='tozero', tickmode="auto")

    del sns_plot
    
    return cutoff_fig

def noisy_rf_difference_density(rf_df_original, rf_df_denoised, name = "Not specified"):
    tmp = pd.DataFrame({ "cutoff" : rf_df_original["RF"] - rf_df_denoised["RF"]})
    sns_plot = sns.kdeplot(data=tmp, x="cutoff")
    del tmp
    
    rf_difference_fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    rf_difference_fig.add_trace(
        go.Histogram(
            x=rf_df_original["RF"] - rf_df_denoised["RF"],
            name='original',
            marker_color='#636EFA',
            nbinsx=50,
            opacity=0.75),
        secondary_y=False)

    rf_difference_fig.update_traces(marker_line_width=1,marker_line_color="black")

    rf_difference_fig.update_layout(
        title= f"Изменение в RF для {name}",
        yaxis_title="Число последовательностей",
        xaxis_title="разность RF",
        template="none",
    )
    
    rf_difference_fig.update_yaxes(title_text="Число последовательностей", rangemode='tozero', secondary_y=False)
    rf_difference_fig.update_yaxes(title_text="Плотность", secondary_y=True, showgrid=False,
                     visible=False, showticklabels=False, rangemode='tozero', tickmode="auto")
    del sns_plot
        
    return rf_difference_fig

def RF_length_KDE(rf_out_df, name = "Not specified"):
    rf_out_df["RF_original_length_log"] = np.log(rf_out_df["RF_original_length"])
    sns_plot = sns.kdeplot(data=rf_out_df, x="RF_original_length_log", hue="RF_bool")
    
    RF_length_fig = go.Figure()
    RF_length_fig.add_trace(go.Scatter(x=sns_plot.get_lines()[0].get_data()[0], 
                             y=sns_plot.get_lines()[0].get_data()[1], 
                             mode="lines", name="Исходный RF выше"))
    RF_length_fig.add_trace(go.Scatter(x=sns_plot.get_lines()[1].get_data()[0], 
                             y=sns_plot.get_lines()[1].get_data()[1], 
                             mode="lines", name="Исходный RF ниже"))
    RF_length_fig.add_trace(go.Scatter(x=sns_plot.get_lines()[2].get_data()[0], 
                             y=sns_plot.get_lines()[2].get_data()[1], 
                             mode="lines", name="Без изменений"))

    RF_length_fig.update_layout(
        title= f"Распределение RF для разных классов в зависимости от длины последовательностей в {name}",
        xaxis_title="Длина последовательности (log)",
        yaxis_title="Плотность",
        template="none",
    )
    
    del sns_plot
    
    return RF_length_fig

def RF_width_KDE(rf_out_df, name = "Not specified"):
    #rf_out_df["RF_original_length_log"] = np.log(rf_out_df["RF_original_length"])
    sns_plot = sns.kdeplot(data=rf_out_df, x="RF_original_width", hue="RF_bool")
    
    RF_width_fig = go.Figure()
    RF_width_fig.add_trace(go.Scatter(x=sns_plot.get_lines()[0].get_data()[0], 
                             y=sns_plot.get_lines()[0].get_data()[1], 
                             mode="lines", name="Исходный RF выше"))
    RF_width_fig.add_trace(go.Scatter(x=sns_plot.get_lines()[1].get_data()[0], 
                             y=sns_plot.get_lines()[1].get_data()[1], 
                             mode="lines", name="Исходный RF ниже"))
    RF_width_fig.add_trace(go.Scatter(x=sns_plot.get_lines()[2].get_data()[0], 
                             y=sns_plot.get_lines()[2].get_data()[1], 
                             mode="lines", name="Без изменений"))

    RF_width_fig.update_layout(
        title= f"Распределение RF для разных классов в зависимости от количества последовательностей в {name}",
        xaxis_title="Количество последовательностей",
        yaxis_title="Плотность",
        template="none",
    )
    
    del sns_plot
        
    return RF_width_fig

def RF_cutoff_KDE(rf_out_df, name):
    #rf_out_df["RF_original_length_log"] = np.log(rf_out_df["RF_original_length"])
    sns_plot = sns.kdeplot(data=rf_out_df, x="cutoff_percent", hue="RF_bool")
    
    RF_cutoff_fig = go.Figure()
    RF_cutoff_fig.add_trace(go.Scatter(x=sns_plot.get_lines()[0].get_data()[0], 
                             y=sns_plot.get_lines()[0].get_data()[1], 
                             mode="lines", name="Исходный RF выше"))
    RF_cutoff_fig.add_trace(go.Scatter(x=sns_plot.get_lines()[1].get_data()[0], 
                             y=sns_plot.get_lines()[1].get_data()[1], 
                             mode="lines", name="Исходный RF ниже"))
    RF_cutoff_fig.add_trace(go.Scatter(x=sns_plot.get_lines()[2].get_data()[0], 
                             y=sns_plot.get_lines()[2].get_data()[1], 
                             mode="lines", name="Без изменений"))

    RF_cutoff_fig.update_layout(
        title= f"Распределение RF для разных классов в зависимости от процента отрезанной длины в {name}",
        xaxis_title="Процент отрезанной длины",
        yaxis_title="Плотность",
        template="none",
    )
    
    del sns_plot
    
    return RF_cutoff_fig

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
    rf_out_df["cutoff_percent"] = (rf_out_df["RF_original_length"] - rf_out_df["RF_denoised_length"])/rf_out_df["RF_original_length"]
    return rf_out_df



