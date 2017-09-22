import matplotlib
matplotlib.use('Agg')
from Swing.util.BoxPlot import BoxPlot
from matplotlib.backends.backend_pdf import PdfPages

import pdb

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import time

def get_df(df, fp, min_lag, max_lag, td_window):
    new_df = df[(df['file_path'] == fp) & (df['min_lag'] == min_lag) & (df['max_lag'] == max_lag) & (df['td_window'] == td_window)]
    return(new_df)

def read_tdr_results(folder_list, folder_str):
    agg_df = pd.DataFrame()
    for input_folder in folder_list:
        for file_path in os.listdir(input_folder):
            if folder_str in file_path:
              df = pd.read_csv(input_folder+file_path,sep='\t', engine='python')
              # check if the columns are misaligned.
              if type(df['permutation_n'].iloc[0]) is str:
                  new_col = df.columns.tolist()
                  new_col.pop(0)
                  new_df = df.iloc[:,0:len(df.iloc[0])-1]
                  new_df.columns = new_col
                  df=new_df

              agg_df = agg_df.append(df)
    return(agg_df)

#input_folder_list = ["/projects/p20519/roller_output/high_sampling/RandomForest/"]  
input_folder_list = ["/projects/p20519/roller_output/gnw/Lasso/"]  
test_statistic = ['aupr', 'auroc']
save_tag = "window_scan"
n_trials = 100

start = time.time()
agg_df = read_tdr_results(input_folder_list, folder_str = "2017-09")
agg_df.to_pickle("RF_window_scan.pkl")
end = time.time()

stat = 'aupr'
network_list = agg_df['file_path'].unique().tolist()

window_sizes = range(1,21)

for network in network_list:
    baseline = get_df(agg_df, network, 0, 0, 21)
    all_available_windows = agg_df[agg_df['file_path'] == network]['td_window'].unique()
    if len(all_available_windows) < len(window_sizes):
        continue
    if len(baseline) == 0:
        continue
    for td_window in window_sizes:
        baseline_mean=baseline[stat].mean()
        if 21-td_window > 2:
            max_lag = 3 
        else:
            max_lag = 21-td_window
        comparisons = get_df(agg_df, network, 1, max_lag, td_window)
        winstat = (comparisons[stat]-baseline_mean)/baseline_mean
        stat_list.append(winstat)

pdb.set_trace()



# group each result by td_window for each network

# final result:
# each dot = each network auroc or aupr/mean of baseline method of each network

completed_files = agg_df['file_path'].unique().tolist()

job_file = pd.read_csv('/home/jjw036/Roller/pipelines/job_params_high_sampling.txt', sep = ' ', header=None)

job_file.columns = ['data_path','data_path', 'input_file', 'iterating_param', 'iterating_style']

mask = job_file['input_file'].str.contains('|'.join(completed_files))

new_jobs = job_file[~mask]

new_jobs.to_csv('/home/jjw036/Roller/pipelines/job_params_RF_window_scan_missing.txt', sep = ' ', header=None, index=False)
print("Added missing jobs")
