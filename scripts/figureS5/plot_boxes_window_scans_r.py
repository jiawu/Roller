import matplotlib
matplotlib.use('Agg')
from Swing.util.BoxPlot import BoxPlot
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

import pdb

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import time
from Swing.util.mplstyle import style1
import seaborn as sns
from palettable.colorbrewer.qualitative import Set1_3

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

input_folder_list = ["/projects/p20519/roller_output/gnw/RandomForest/"]  
test_statistic = ['aupr', 'auroc']
save_tag = "window_scan"
n_trials = 100

start = time.time()
agg_df = read_tdr_results(input_folder_list, folder_str = "2017-09")
#agg_df.to_pickle("Dionesus_window_scan.pkl")
#agg_df = pd.read_pickle("Dionesus_window_scan.pkl")
end = time.time()

stat = 'aupr'
network_list = agg_df['file_path'].unique().tolist()

window_sizes = range(1,22)
outer_list = []
overall_df = pd.DataFrame()

for td_window in window_sizes:
    inner_list = []
    for network in network_list:
        baseline = get_df(agg_df, network, 0, 0, 21)
        if len(baseline) == 0:
            continue
        if 21-td_window > 2:
            max_lag = 3 
        else:
            max_lag = 21-td_window
        if (td_window == 21):
            min_lag = 0
            max_lag = 0
        else:
            min_lag = 1
        
        comparisons = get_df(agg_df, network, min_lag, max_lag, td_window)
        if len(comparisons) == 0:
            continue
        stat = 'aupr'
        baseline_mean=baseline[stat].mean()
        comparisons['percent_{}'.format(stat)] = ((comparisons[stat]-baseline_mean)/baseline_mean)*100
        stat = 'auroc'
        baseline_mean=baseline[stat].mean()
        comparisons['percent_{}'.format(stat)] = ((comparisons[stat]-baseline_mean)/baseline_mean)*100
        overall_df = overall_df.append(comparisons.iloc[0:50,:], ignore_index = True)
    outer_list.append(inner_list)

stat = 'percent_aupr'

colors = []
for w in range(1, 21):
    test_data = overall_df[overall_df.td_window == w]
    baseline = overall_df[overall_df.td_window == 21]
    baseline_mean = baseline[stat].mean()
    diff = np.mean(test_data[stat])-baseline_mean
    if stats.ttest_ind(test_data[stat], baseline[stat])[1] < 0.05:
        if diff > 0:
            colors.append(Set1_3.mpl_colors[0])
        else:
            colors.append(Set1_3.mpl_colors[1])
    else:
        colors.append('grey')

fig, ax = plt.subplots(figsize=(11,7))
sns.boxplot(ax = ax, data = overall_df, x = 'td_window', y = 'percent_aupr', palette=colors)
xlabs = ax.get_xticks()
ax.set_xticklabels(['{:d}'.format(x+1) for x in xlabs])
ax.set_ylabel('Percent Difference AUPR')
ax.set_xlabel('Window Size')
fig.savefig('RandomForest_10_AUPR_window_scan.png')

fig, ax = plt.subplots(figsize=(11,7))
sns.boxplot(ax = ax, data = overall_df, x = 'td_window', y = 'percent_auroc', palette=colors)
xlabs = ax.get_xticks()
ax.set_xticklabels(['{:d}'.format(x+1) for x in xlabs])
ax.set_ylabel('Percent Difference AUROC')
ax.set_xlabel('Window Size')
fig.savefig('RandomForest_10_AUROC_window_scan.png')

