import matplotlib
matplotlib.use('Agg')
from Roller.util.BoxPlot import BoxPlot
import pdb

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

"""
Script that loads data from a dataframe and generates boxplots

"""

output_path = "/home/jjw036/"
input_folder = "/projects/p20519/roller_output/stability_analysis/RandomForest/"
test_statistic = 'aupr'
save_tag = "Insilico_aupr_n_perms_stability"
n_trials = 200

agg_df = pd.DataFrame()
for file_path in os.listdir(input_folder):    
    df = pd.read_csv(input_folder+file_path,sep=',|\t', engine='python')
    agg_df = agg_df.append(df)

datasets = ['insilico_size10_1','insilico_size10_2','insilico_size10_3','insilico_size10_4','insilico_size10_5']

label_list = []
auroc_list = []

#stability analysis
for dataset in datasets:
    current_df = agg_df[agg_df['file_path'].str.contains(dataset)]
    test_stat_10 = current_df[(current_df['n_trees'] ==500) & (current_df['permutation_n']==10)]
    test_stat_100 = current_df[(current_df['n_trees'] ==500) & (current_df['permutation_n']==100)]
    test_stat_500 = current_df[(current_df['n_trees'] ==500) & (current_df['permutation_n']==500)]
    test_stat_1000 = current_df[(current_df['n_trees'] ==500) & (current_df['permutation_n']==1000)]
    pdb.set_trace()
    auroc_list.append(test_stat_10[test_statistic][0:n_trials].tolist())
    auroc_list.append(test_stat_100[test_statistic][0:n_trials].tolist())
    auroc_list.append(test_stat_500[test_statistic][0:n_trials].tolist())
    auroc_list.append(test_stat_1000[test_statistic][0:n_trials].tolist())
    
    label_list.append("n_perm = 10")
    label_list.append("n_perm = 100")
    label_list.append("n_perm = 500")
    label_list.append("n_perm = 1000")


bp_data = auroc_list
bp = BoxPlot()
bp.plot_box(bp_data, label_list)
#auroc_1 = df['auroc'].values
#auroc_2 = df['auroc'].values

#bp_data = [auroc_1,auroc_2]

#bp = BoxPlot()

#bp.plot_box(bp_data, ['n_trees = 10', 'n_trees = 20'])

bp.add_formatting()        

bp.save_plot(output_path, save_tag)

#grouped.get_group((2,2)).mean()['aupr']


