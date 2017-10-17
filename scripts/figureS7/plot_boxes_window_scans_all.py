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

def get_df(df, fp, min_lag, max_lag, td_window, inftype = "RandomForest"):
    new_df = df[(df['file_path'] == fp) & (df['min_lag'] == min_lag) & (df['max_lag'] == max_lag) & (df['td_window'] == td_window) & (df['InfType'] == inftype)]
    return(new_df)

def load_data():
    input_folder_list = ["/projects/p20519/roller_output/gnw/RandomForest/"]  
    agg_df_RF = read_tdr_results(input_folder_list, folder_str = "2017-09")
    agg_df_RF['InfType'] = 'RandomForest'

    input_folder_list = ["/projects/p20519/roller_output/gnw/Dionesus/"]  
    agg_df_P = read_tdr_results(input_folder_list, folder_str = "2017-09")
    agg_df_P['InfType'] = 'PLSR'

    input_folder_list = ["/projects/p20519/roller_output/gnw/Lasso/"]  
    agg_df_L = read_tdr_results(input_folder_list, folder_str = "2017-09")
    agg_df_L['InfType'] = 'Lasso'

    all_dfs = [agg_df_RF, agg_df_P, agg_df_L]
    merged_df = pd.concat(all_dfs)
    return(merged_df)
    

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

def get_inf_df(network_1, inf_type):
    RFnet1 = network_1[network_1['InfType'] == inf_type]
    RFn1 = RFnet1.groupby('td_window').mean()
    return(RFn1)



def get_comparisons(merged_df, inftypes, window_sizes, network_list):
    overall_df = pd.DataFrame()
    network_1_df = pd.DataFrame()
    for inftype in inftypes:
        for td_window in window_sizes:
            for network in network_list:
                baseline = get_df(merged_df, network, 0, 0, 21, inftype = inftype)
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
                
                comparisons = get_df(merged_df, network, min_lag, max_lag, td_window, inftype = inftype)
                if len(comparisons) == 0:
                    continue
                
                # for each statistic, get the percent difference to the baseline comparison.
                stat = 'aupr'
                baseline_mean=baseline[stat].mean()
                comparisons['percent_{}'.format(stat)] = ((comparisons[stat]-baseline_mean)/baseline_mean)*100
                stat = 'auroc'
                baseline_mean=baseline[stat].mean()
                comparisons['percent_{}'.format(stat)] = ((comparisons[stat]-baseline_mean)/baseline_mean)*100
                overall_df = overall_df.append(comparisons.iloc[0:20,:], ignore_index = True)
                if network == network_list[6]:
                    network_1_df = network_1_df.append(comparisons.iloc[0:20,:], ignore_index = True)
                print(comparisons,len(comparisons))
    return(overall_df, network_1_df)

test_statistic = ['aupr', 'auroc']
save_tag = "window_scan"
n_trials = 100

#merged_df = load_data()
#merged_df.to_pickle("merged_window_scan.pkl")

#merged_df = pd.read_pickle("merged_window_scan.pkl")

#network_list = merged_df['file_path'].unique().tolist()
#window_sizes = range(2,22)
#inftypes = ['RandomForest', 'Lasso', 'PLSR']

#overall_df, network_1 = get_comparisons(merged_df, inftypes, window_sizes, network_list)

#overall_df.to_pickle("merged_window_scan_comparisons.pkl")
#network_1.to_pickle("merged_window_scan_comparisons_network1.pkl")
overall_df = pd.read_pickle("merged_window_scan_comparisons.pkl")
network_1 = pd.read_pickle("merged_window_scan_comparisons_network1.pkl")

stat = 'percent_aupr'

meanpointprops = dict(marker = 'o', markeredgecolor='black', markerfacecolor='white', markeredgewidth = 1, markersize = 10)
medianprops = dict(color='black', linewidth=1.5)
g = sns.FacetGrid(overall_df, col = "InfType", size = 12, aspect=0.6)
g = g.map(sns.boxplot, 'td_window', 'percent_aupr', showfliers = False, color = 'cornflowerblue', showmeans=True, meanprops = meanpointprops, medianprops=medianprops)
g.set(ylim=(-60, 200))
g.set(xlabel="Window Size")
g.fig.get_axes()[0].set_ylabel("% change AUPR")
g.fig.get_axes()[0].axhline(y = 0, c = "darkgrey")
g.fig.get_axes()[1].axhline(y = 0, c = "darkgrey")
g.fig.get_axes()[2].axhline(y = 0, c = "darkgrey")

## Plot the AUPR for 1 graph, dionesus
RFn1 = get_inf_df(network_1, inf_type = 'RandomForest')
net1color = 'firebrick'
g.fig.get_axes()[0].plot(RFn1[stat].values, marker='.', color = net1color)

Ln1 = get_inf_df(network_1, inf_type = 'Lasso')
g.fig.get_axes()[1].plot(Ln1[stat].values, marker='.', color = net1color)

Pn1 = get_inf_df(network_1, inf_type = 'PLSR')
g.fig.get_axes()[2].plot(Pn1[stat].values, marker='.', color = net1color)

for axes in g.axes.flat:
    axes.set_xticklabels(['{:d}'.format(x+2) for x in axes.get_xticks()])
g.savefig('combined_10_AUPR_window_scan.pdf')

## Plot AUROC

stat = 'percent_auroc'

meanpointprops = dict(marker = 'o', markeredgecolor='black', markerfacecolor='white', markeredgewidth = 1, markersize = 10)
medianprops = dict(color='black', linewidth=1.5)
g = sns.FacetGrid(overall_df, col = "InfType", size = 12, aspect=0.6)
g = g.map(sns.boxplot, 'td_window', 'percent_auroc', showfliers = False, color = 'darkorange', showmeans=True, meanprops = meanpointprops, medianprops=medianprops)
g.set(ylim=(-50, 60))
g.set(xlabel="Window Size")
g.fig.get_axes()[0].set_ylabel("% change AUROC")
g.fig.get_axes()[0].axhline(y = 0, c = "darkgrey")
g.fig.get_axes()[1].axhline(y = 0, c = "darkgrey")
g.fig.get_axes()[2].axhline(y = 0, c = "darkgrey")

## Plot the AUROC for 1 graph, all inf methods
RFn1 = get_inf_df(network_1, inf_type = 'RandomForest')
net1color = 'firebrick'
g.fig.get_axes()[0].plot(RFn1[stat].values, marker='.', color = net1color)

Ln1 = get_inf_df(network_1, inf_type = 'Lasso')
g.fig.get_axes()[1].plot(Ln1[stat].values, marker='.', color = net1color)

Pn1 = get_inf_df(network_1, inf_type = 'PLSR')
g.fig.get_axes()[2].plot(Pn1[stat].values, marker='.', color = net1color)

for axes in g.axes.flat:
    axes.set_xticklabels(['{:d}'.format(x+2) for x in axes.get_xticks()])

g.savefig('combined_10_AUROC_window_scan.pdf')

