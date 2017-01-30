import sys
from datetime import datetime
import numpy as np

sys.path.append("../pipelines")
import Pipelines as tdw
import Swing.util.lag_identification as lag_id
from Swing.util.Evaluator import Evaluator

data_folder = "../data/invitro/Dionesus"
file_path = "../data/invitro/omranian_parsed_timeseries.tsv"
#file_path = "../data/invitro/cantone_switchon_interpolated_timeseries.tsv"
#file_path = "../data/gnw_insilico/network_data/Ecoli/Ecoli-9_timeseries.tsv"
current_time = datetime.time
run_params = {'data_folder': data_folder,
              'file_path':file_path,
              'td_window':17,
              'min_lag':1,
              'max_lag':1,
              'n_trees':100,
              'permutation_n':10,
              'lag_method':'mean_mean',
              'calc_mse':False,
              'bootstrap_n':1000,
              'n_trials':1,
              'run_time':current_time,
              'sort_by':'rank',
              'iterating_param':'td_window',
              'filter_noisy': False
              }

window_sizes = list(range(3,6))
roc_list = []
pr_list = []

for window_size in window_sizes:
    print("window_size:",window_size)
    run_params['td_window'] = window_size
    if window_size == 5:
        run_params['min_lag'] =0
        run_params['max_lag'] = 0
    roc,pr, tdr = tdw.get_td_stats_custom(**run_params)
    roc_list.append(roc)
    pr_list.append(pr)
