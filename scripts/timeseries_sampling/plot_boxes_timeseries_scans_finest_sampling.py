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

def get_df(df, sampling_rate, network, min_lag, max_lag, td_window, inftype = "RandomForest"):
    new_df = df[(df['network_number'] == network) & (df['sampling_rate'] == sampling_rate) & (df['min_lag'] == min_lag) & (df['max_lag'] == max_lag) & (df['td_window'] == float(td_window)) & (df['InfType'] == inftype)]
    return(new_df)

def load_data():
    input_folder_list = ["/projects/p20519/roller_output/high_sampling/RandomForest/"]  
    agg_df_RF = read_tdr_results(input_folder_list, folder_str = "2017-09")
    agg_df_RF['InfType'] = 'RandomForest'

    input_folder_list = ["/projects/p20519/roller_output/high_sampling/Dionesus/"]  
    agg_df_P = read_tdr_results(input_folder_list, folder_str = "2017-09")
    agg_df_P['InfType'] = 'PLSR'

    input_folder_list = ["/projects/p20519/roller_output/high_sampling/Lasso/"]  
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
              try:
                  df = pd.read_csv(input_folder+file_path,sep='\t', engine='python')
              except pd.io.common.EmptyDataError:
                  continue
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
    RFn1 = RFnet1.groupby('sampling_rate').mean()
    return(RFn1)


def get_comparisons(merged_df, inftypes, sampling_rates, network_list):
    # 10 - 101, 67
    # 30 - 34, 22
    # 50 - 21, 14
    # 100 - 11, 7
    # 200 - 6, 4
    # 333 - 4, 3
    # 500 - 3, 2

    max_window_dict = { '10':101,
                        '30':34,
                        '50':21,
                        '100':11,
                        '200':6,
                        '333':4,
                        '500':3}
    td_window_dict = {  '10':67,
                        '30':22,
                        '50':14,
                        '100':7,
                        '200':4,
                        '333':2,
                        '500':2}
    all_networks = merged_df['network_number'].dropna().unique()
    all_sampling_rates = merged_df['sampling_rate'].dropna().unique()
    overall_df = pd.DataFrame()
    network_1_df = pd.DataFrame()
    for inftype in inftypes:
        for sampling_rate in all_sampling_rates:
            for network in all_networks:
                baseline = get_df(merged_df, sampling_rate, network, 0, 0, max_window_dict[sampling_rate], inftype = inftype)
                if len(baseline) == 0:
                    continue
                td_window = td_window_dict[sampling_rate]
                min_lag = 1
                max_lag = 3
                max_window = max_window_dict[sampling_rate]
                if max_window - td_window < 3:
                    max_lag = max_window - td_window
                    if min_lag > max_lag:
                        min_lag = max_lag
                comparisons = get_df(merged_df, sampling_rate, network, min_lag, max_lag, td_window, inftype = inftype)
                if len(comparisons) == 0:
                    continue
                
                # for each statistic, get the percent difference to the baseline comparison.
                stat = 'aupr'
                baseline_mean=baseline[stat].mean()
                comparisons['percent_{}'.format(stat)] = ((comparisons[stat]-baseline_mean)/baseline_mean)*100
                stat = 'auroc'
                baseline_mean=baseline[stat].mean()
                comparisons['percent_{}'.format(stat)] = ((comparisons[stat]-baseline_mean)/baseline_mean)*100
                overall_df = overall_df.append(comparisons.iloc[0:100,:], ignore_index = True)
                if network == all_networks[6]:
                    network_1_df = network_1_df.append(comparisons.iloc[0:100,:], ignore_index = True)
                print(sampling_rate, network, len(baseline),len(comparisons), inftype)
    return(overall_df, network_1_df)

test_statistic = ['aupr', 'auroc']
save_tag = "window_scan"
n_trials = 100

#merged_df = load_data()
#merged_df.to_pickle("merged_sampling_scan.pkl")

#merged_df = pd.read_pickle("merged_sampling_scan.pkl")

#network_list = merged_df['file_path'].unique().tolist()
#network_list = [x for x in network_list if 'even' not in str(x)]

#merged_df['network_number'] = merged_df['file_path'].str.extract('((Ecoli|Yeast)-[0-9]{1,2})_[\d]+_timeseries')[0]
#merged_df['sampling_rate'] = merged_df['file_path'].str.extract('-[0-9]{1,2}_(\d+)_timeseries')

# baseline is the full size window. The full size windows are:
# 
# parse the dataframe:
# sampling period
# base_network

# filter data

window_sizes = range(2,22)
inftypes = ['RandomForest', 'Lasso', 'PLSR']
sampling_rates = [10]
overall_df, network_1 = get_comparisons(merged_df, inftypes, sampling_rates, network_list)

overall_df.to_pickle("merged_sampling_comparisons_finest.pkl")
network_1.to_pickle("merged_sampling_comparisons_network1_finest.pkl")
#overall_df = pd.read_pickle("merged_sampling_comparisons_finest.pkl")
#network_1 = pd.read_pickle("merged_sampling_comparisons_network1_finest.pkl")

stat = 'percent_aupr'
pdb.set_trace()
meanpointprops = dict(marker = 'o', markeredgecolor='black', markerfacecolor='white', markeredgewidth = 1, markersize = 10)
medianprops = dict(color='black', linewidth=1.5)
g = sns.FacetGrid(overall_df, col = "InfType", size = 12, aspect=0.6)
g = g.map(sns.boxplot, 'sampling_rate', 'percent_aupr', showfliers = False, order = ['10','30','50','100','200','333','500'],color = 'cornflowerblue', showmeans=True, meanprops = meanpointprops, medianprops=medianprops)
g.set(ylim=(-60, 200))
g.set(xlabel="Sampling Rate")
g.fig.get_axes()[0].set_ylabel("% change AUPR")
g.fig.get_axes()[0].axhline(y = 0, c = "darkgrey")
g.fig.get_axes()[1].axhline(y = 0, c = "darkgrey")
g.fig.get_axes()[2].axhline(y = 0, c = "darkgrey")

## Plot the AUPR for 1 graph, dionesus
index = ['10', '30', '50', '100', '200','333', '500']
RFn1 = get_inf_df(network_1, inf_type = 'RandomForest')
net1color = 'firebrick'
g.fig.get_axes()[0].plot(RFn1.loc[index, stat], marker='.', color = net1color)

Ln1 = get_inf_df(network_1, inf_type = 'Lasso')
g.fig.get_axes()[1].plot(Ln1.loc[index, stat], marker='.', color = net1color)

Pn1 = get_inf_df(network_1, inf_type = 'PLSR')
g.fig.get_axes()[2].plot(Pn1.loc[index,stat], marker='.', color = net1color)

#for axes in g.axes.flat:
#    axes.set_xticklabels(['{:d}'.format(x) for x in axes.get_xticks()])
g.savefig('combined_10_AUPR_sampling_heat_map.pdf')

pdb.set_trace()
## Plot AUROC

stat = 'percent_auroc'

meanpointprops = dict(marker = 'o', markeredgecolor='black', markerfacecolor='white', markeredgewidth = 1, markersize = 10)
medianprops = dict(color='black', linewidth=1.5)
g = sns.FacetGrid(overall_df, col = "InfType", size = 12, aspect=0.6)
g = g.map(sns.boxplot, 'sampling_rate', 'percent_auroc', order=['10','30','50','100','200','333','500'], showfliers = False, color = 'darkorange', showmeans=True, meanprops = meanpointprops, medianprops=medianprops)
g.set(ylim=(-50, 60))
g.set(xlabel="Sampling Rate")
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

#for axes in g.axes.flat:
#    axes.set_xticklabels(['{:d}'.format(x) for x in axes.get_xticks()])

g.savefig('combined_10_AUROC_sampling_heatmap.pdf')

