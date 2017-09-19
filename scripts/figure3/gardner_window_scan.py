# load packages
import pandas as pd
import statsmodels.tsa.stattools as stats
import statsmodels.graphics.tsaplots as sg
import matplotlib.pyplot as plt
import matplotlib
import sys
from datetime import datetime
import numpy as np
import warnings
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
run_params = {'data_folder': data_folder,
              'file_path':file_path,
              'td_window':10,
              'min_lag':0,
              'max_lag':1,
              'n_trees':1000,
              'permutation_n':None,
              'lag_method':'mean_mean',
              'calc_mse':False,
              'bootstrap_n':1000,
              'n_trials':1,
              'run_time':current_time,
              'sort_by':'rank',
              'window_type': 'RandomForest'
              }

preexisting = False
if not preexisting:
    roc_list =[]
    pr_list = []
    rankings = []
    for ii in range(50):
        print("Run: ", str(ii))
        roc, pr, tdr = tdw.get_td_stats(**run_params)
        roc_list.append(roc)
        pr_list.append(pr)
        # rankings.append(edge_list)