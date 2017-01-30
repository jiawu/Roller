import matplotlib
matplotlib.use('Agg')
from Swing.util.BoxPlot import BoxPlot
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

input_folder_list = ["/projects/p20519/roller_output/stability_analysis/RandomForest/","/projects/p20519/roller_output/stability_analysis/Lasso/","/projects/p20519/roller_output/stability_analysis/Dionesus/"]  
test_statistic = 'aupr'
save_tag = "Insilico_methods_comparison_aupr"
n_trials = 50

agg_df = pd.DataFrame()
for input_folder in input_folder_list:
    for file_path in os.listdir(input_folder):    
        df = pd.read_csv(input_folder+file_path,sep=',|\t', engine='python')
        agg_df = agg_df.append(df)

datasets = ['insilico_size10_1','insilico_size10_2','insilico_size10_3','insilico_size10_4','insilico_size10_5']

label_list = []
auroc_list = []

#stability analysis
for dataset in datasets[0:1]:
    current_df = agg_df[agg_df['file_path'].str.contains(dataset)]
    test_stat_D_21 = current_df[(current_df['data_folder'].str.contains('Dionesus')) & (current_df['td_window'] == 21)]
    test_stat_L_21 = current_df[(current_df['data_folder'].str.contains('Lasso')) & (current_df['td_window'] == 21)]
    test_stat_RF_21 = current_df[(current_df['data_folder'].str.contains('RandomForest')) & (current_df['td_window'] == 21)]
    pdb.set_trace()
    
    test_stat_D_10 = current_df[(current_df['data_folder'].str.contains('Dionesus')) & (current_df['td_window'] == 10)]
    test_stat_L_10 = current_df[(current_df['data_folder'].str.contains('Lasso')) & (current_df['td_window'] == 10)]
    test_stat_RF_10 = current_df[(current_df['data_folder'].str.contains('RandomForest')) & (current_df['td_window'] == 10)]

    test_stat_D_15 = current_df[(current_df['data_folder'].str.contains('Dionesus')) & (current_df['td_window'] == 15)]
    test_stat_L_15 = current_df[(current_df['data_folder'].str.contains('Lasso')) & (current_df['td_window'] == 15)]
    test_stat_RF_15 = current_df[(current_df['data_folder'].str.contains('RandomForest')) & (current_df['td_window'] == 15)]


    auroc_list.append(test_stat_D_21[test_statistic][0:n_trials].tolist())
    auroc_list.append(test_stat_L_21[test_statistic][0:n_trials].tolist())
    auroc_list.append(test_stat_RF_21[test_statistic][0:n_trials].tolist())
    
    auroc_list.append(test_stat_D_10[test_statistic][0:n_trials].tolist())
    auroc_list.append(test_stat_L_10[test_statistic][0:n_trials].tolist())
    auroc_list.append(test_stat_RF_10[test_statistic][0:n_trials].tolist())
    
    auroc_list.append(test_stat_D_15[test_statistic][0:n_trials].tolist())
    auroc_list.append(test_stat_L_15[test_statistic][0:n_trials].tolist())
    auroc_list.append(test_stat_RF_15[test_statistic][0:n_trials].tolist())
    
    label_list.append("Dionesus")
    label_list.append("Lasso")
    label_list.append("RandomForest")
    label_list.append("td-Dionesus, w=10")
    label_list.append("td-Lasso, w=10")
    label_list.append("td-RandomForest, w=10")
    label_list.append("td-Dionesus, w=15")
    label_list.append("td-Lasso, w=15")
    label_list.append("td-RandomForest, w=15")


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


