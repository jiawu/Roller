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

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def get_df(df, fp, min_lag, max_lag, td_window):
    new_df = df[(df['file_path'] == fp) & (df['min_lag'] == min_lag) & (df['max_lag'] == max_lag) & (df['td_window'] == td_window)]
    return(new_df)

def get_df_mm(df, fp, min_lag, max_lag):
    new_df = df[(df['file_path'] == fp) & (df['min_lag'] == min_lag) & (df['max_lag'] == max_lag)]
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
input_folder_list = ["/projects/p20519/roller_output/gnw/Dionesus/"]  
test_statistic = ['aupr', 'auroc']
save_tag = "window_scan"
n_trials = 100

start = time.time()
agg_df = read_tdr_results(input_folder_list, folder_str = "2017-09")
#agg_df.to_pickle("RF_window_scan.pkl")
end = time.time()

stat = 'aupr'
network_list = agg_df['file_path'].unique().tolist()

window_sizes = range(1,21)
organism = "aggregate"
inf_method = "Dionesus"

outer_list = []

min_maxes = [(0,0), (0,1), (0,2), (0,3), (0,4), (1,1), (1,2), (1,3), (1,4), (2,2), (2,3), (2,4), (3,3), (3,4), (4,4)]

for mm in min_maxes:
    inner_list = []
    for network in network_list:
        #if organism not in network:
        #    continue
        baseline = get_df(agg_df, network, 0, 0, 21)
        if len(baseline) == 0:
            continue
        baseline_mean=baseline[stat].mean()
        
        td_window = 15
        max_lag = mm[1]
        if td_window + max_lag > 21:
            td_window = 21-max_lag
        comparisons = get_df(agg_df, network, min_lag = mm[0], max_lag = mm[1], td_window=td_window)

        if len(comparisons) == 0:
            continue
        
        winstat = ((comparisons[stat]-baseline_mean)/baseline_mean)*100
        inner_list.append(winstat.mean())
    outer_list.append(inner_list)
fig, ax = plt.subplots(figsize=(8,5))
violin = ax.violinplot(outer_list, points = 100, widths = 0.4, showmeans=False, showextrema=False, showmedians=False)

color = "#E49E22"
#color = "#D36228"
#color = "#6CB1D2"
for element in violin['bodies']:
    element.set(color=color)
    element.set(alpha=1)
    element.set(edgecolor='white')
#### Box plot within violinplots

quartile1 = []
medians = []
quartile3 =[]

qmq = []
#qmq.append(tuple(np.percentile(tt1.values, [25,50,75])))
for ws in outer_list:
    qmq.append(tuple(np.percentile(ws, [25,50,75])))

quartile1 = [x[0] for x in qmq]
medians = [x[1] for x in qmq]
quartile3 = [x[2] for x in qmq]

whiskers = np.array([adjacent_values(sorted_array, q1, q3) for sorted_array, q1, q3 in zip(outer_list, quartile1, quartile3)])
whiskersMin, whiskersMax = whiskers[:, 0], whiskers[:, 1]

inds = np.arange(1, len(medians) + 1)
ax.scatter(inds, medians, marker='o', color='white', s=15, zorder=3)
#ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
#ax.vlines(inds, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)
ax.axhline(y=0, linestyle='-', color='darkgrey')

ax.set_ylabel('Percentage Increase of {} (PLSR)'.format(stat.upper()))
ax.set_xlabel('Min/Max')

ax.set_ylim([-50,150])
ax.set_xlim([ 0, len(min_maxes)+1])
ax.set_yticks([x for x in range(-50,150,25)])
ax.set_xticks([x for x in range(1, len(min_maxes)+1)])
ax.set_xticklabels([str(x) for x in min_maxes])

fp = "{}_{}_{}.pdf".format(organism,inf_method, "min_max")
fig.tight_layout()
plt.savefig(fp)


