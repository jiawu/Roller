# load packages
import pandas as pd
import statsmodels.tsa.stattools as stats
import statsmodels.graphics.tsaplots as sg
import matplotlib.pyplot as plt
import matplotlib
import itertools as it
import sys
from datetime import datetime
import numpy as np
import warnings
import json
import time

warnings.filterwarnings("ignore")

import networkx as nx
from nxpd import draw
from nxpd import nxpdParams

nxpdParams['show'] = 'ipynb'

sys.path.append("../../pipelines")
import Pipelines as tdw

data_folder = "../../data/invitro/"

output_path = "../../data/invitro/"

current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

file_path = "../../data/invitro/gardner_timeseries.tsv"

n_trials = 50
run_params = {'data_folder': data_folder,
              'file_path': file_path,
              'td_window': 7,
              'min_lag': 0,
              'max_lag': 1,
              'n_trees': 1000,
              'permutation_n': 0,
              'lag_method': 'mean_mean',
              'calc_mse': False,
              'bootstrap_n': 1000,
              'n_trilas50': n_trials,
              'run_time': current_time,
              'sort_by': 'rank',
              'window_type': 'RandomForest'
              }

overall_df = pd.DataFrame()
for n in range(n_trials):
    print('Trial {} \n=====\n====='.format(n))
    roc, pr, tdr, edge_list = tdw.get_td_stats(**run_params)
    run_params['roc'] = roc
    run_params['pr'] = pr
    run_params['edge_list'] = edge_list
    run_result = pd.Series(run_params)
    overall_df = overall_df.append(run_result, ignore_index=True)
overall_df.to_pickle('gardner_RF_w7_results.pkl')
sys.exit()


# test = str(list(it.combinations_with_replacement(range(0,5), 2))).strip('[]')
# test = test.replace('(', "").replace(')', "").replace(', ', '_')
# print(test)
# sys.exit()

roc_list = []
pr_list = []
for w in range(2,15):
    run_params['td_window'] = w
    if w == 14:
        run_params['max_lag']=0
    roc, pr, tdr, edge_list = tdw.get_td_stats(**run_params)
    roc_list.append(roc)
    pr_list.append(pr)

plt.plot(range(2,14), roc_list[:-1], '.-')
plt.plot([2,15], [roc_list[-1], roc_list[-1]], '-k', zorder=0)
plt.ylim([0, 1])
plt.xlim([2, 13])
plt.ylabel('AUROC')
plt.xlabel('Window Size')

plt.figure()
plt.plot(range(2,14), pr_list[:-1], '.-')
plt.plot([2,15], [pr_list[-1], pr_list[-1]], '-k', zorder=0)
plt.ylim([0, 1])
plt.xlim([2, 13])
plt.ylabel('AUPR')
plt.xlabel('Window Size')
plt.show()
sys.exit()

methods = ['RandomForest']
timepoints = 14
windowsize = range(2, timepoints+1)
n_windows = timepoints-np.array(windowsize)+1
# min_max = list(it.combinations_with_replacement(range(0,timepoints-1), 2))
min_max = list(it.combinations_with_replacement(range(0, 5), 2))

# preexisting = False
# if not preexisting:
#     roc_list = []
#     pr_list = []
#     rankings = []
#     param_list = []
#     for ii in range(1):
#         print("Run: ", str(ii))
#         roc, pr, tdr, edge_list = tdw.get_td_stats(**run_params)
#         roc_list.append(roc)
#         pr_list.append(pr)
#         rankings.append(edge_list)
#         param_list.append(run_params)

results = []


i = 1
num_test = 5
for m in methods:
    cur_params = run_params.copy()
    cur_params['window_type'] = m
    for ws in windowsize:
        cur_params['td_window'] = ws
        for min_lag, max_lag in min_max:
            cur_params['min_lag'] = min_lag
            cur_params['max_lag'] = max_lag

            try:
                roc, pr, tdr, edge_list = tdw.get_td_stats(**cur_params)
            except ValueError:
                roc = 0
                pr = 0
                edge_list = []

            results.append([m, ws, min_lag, max_lag, roc, pr, edge_list])
            print([m, ws, min_lag, max_lag, roc, pr, edge_list])
            i += 1

results_df = pd.DataFrame(results, columns=['method', 'window_size', 'min_lag', 'max_lag', 'roc', 'pr', 'edge_list'])
results_df.to_pickle('dionesus_scan.pkl')