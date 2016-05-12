# load packages
import pandas as pd
import statsmodels.tsa.stattools as stats
import statsmodels.graphics.tsaplots as sg
import matplotlib.pyplot as plt
import matplotlib
import sys
from datetime import datetime
import numpy as np

import networkx as nx

sys.path.append("../pipelines")
import Pipelines as tdw
data_folder = "/projects/p20519/roller_output/large_networks/Lasso/insilico_size100_5/"

output_path = "/home/jjw036/Roller/insilico_size10_1"

current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')

file_path = "../data/gnw_insilico/network_data/Yeast100/Yeast100-5_timeseries.tsv"
alpha_list = np.arange(0.1,1, 0.05)
roc_list = []
pr_list = []

for alpha in alpha_list:
    run_params = {'data_folder': data_folder,
                  'file_path':file_path,
                  'td_window':10,
                  'min_lag':1,
                  'max_lag':3,
                  'n_trees':10,
                  'permutation_n':10,
                  'lag_method':'mean_mean',
                  'calc_mse':False,
                  'bootstrap_n':1000,
                  'n_trials':1,
                  'run_time':current_time,
                  'sort_by':'rank',
                  'iterating_param':'alpha',
                  'alpha': alpha,
                  'filter_noisy': False
                  }

    roc,pr, tdr = tdw.get_td_stats(**run_params)
    roc_list.append(roc)
    pr_list.append(pr)

plt.plot(alpha_list,roc_list) 
plt.savefig('lasso_alpha_scan.png')
