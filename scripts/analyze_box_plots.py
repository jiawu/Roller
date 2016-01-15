import matplotlib
matplotlib.use('Agg')
from Swing.util.BoxPlot import BoxPlot
import pdb

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

def main(dataset_number,test_statistic): 
    #stability analysis
    datasets = ['insilico_size10_1','insilico_size10_2','insilico_size10_3','insilico_size10_4','insilico_size10_5']
    auroc_list = []
    label_list = []
    for dataset in datasets[dataset_number-1:dataset_number]:
        current_df = agg_df[agg_df['file_path'].str.contains(dataset)]
        #current_df = agg_df
        test_stat_D_21 = current_df[(current_df['data_folder'].str.contains('Dionesus')) & (current_df['td_window'] == 21)]
        test_stat_L_21 = current_df[(current_df['data_folder'].str.contains('Lasso')) & (current_df['td_window'] == 21)]
        test_stat_RF_21 = current_df[(current_df['data_folder'].str.contains('RandomForest')) & (current_df['td_window'] == 21)]
        

        test_stat_D_15 = current_df[(current_df['data_folder'].str.contains('Dionesus')) & (current_df['td_window'] == 15)]
        test_stat_L_15 = current_df[(current_df['data_folder'].str.contains('Lasso')) & (current_df['td_window'] == 15)]
        test_stat_RF_15 = current_df[(current_df['data_folder'].str.contains('RandomForest')) & (current_df['td_window'] == 15)]


        auroc_list.append(test_stat_D_21[test_statistic][0:n_trials].tolist())
        auroc_list.append(test_stat_D_15[test_statistic][0:n_trials].tolist())
        
        auroc_list.append(test_stat_L_21[test_statistic][0:n_trials].tolist())
        auroc_list.append(test_stat_L_15[test_statistic][0:n_trials].tolist())
        
        auroc_list.append(test_stat_RF_21[test_statistic][0:n_trials].tolist())
        auroc_list.append(test_stat_RF_15[test_statistic][0:n_trials].tolist())
        
        
        
        label_list.append("Dionesus")
        label_list.append("td-Dionesus")
        
        label_list.append("Lasso")
        label_list.append("td-Lasso")
        
        label_list.append("RandomForest")
        label_list.append("td-RandomForest")


    bp_data = auroc_list
    bp = BoxPlot()
    bp.plot_box(bp_data, label_list)
    #auroc_1 = df['auroc'].values
    #auroc_2 = df['auroc'].values

    #bp_data = [auroc_1,auroc_2]

    #bp = BoxPlot()

    #bp.plot_box(bp_data, ['n_trees = 10', 'n_trees = 20'])

    bp.add_formatting(test_statistic, "Network "+ str(dataset_number))        

    bp.save_plot(output_path, save_tag+"_"+test_statistic+"_network"+str(dataset_number))


output_path = "/home/jjw036/"

input_folder_list = ["/projects/p20519/roller_output/stability_analysis/RandomForest/","/projects/p20519/roller_output/stability_analysis/Lasso/","/projects/p20519/roller_output/stability_analysis/Dionesus/"]  
test_statistic = 'auroc'
save_tag = "Insilico_methods_comparison"
n_trials = 50

agg_df = pd.DataFrame()
for input_folder in input_folder_list:
    for file_path in os.listdir(input_folder):    
        df = pd.read_csv(input_folder+file_path,sep=',|\t', engine='python')
        agg_df = agg_df.append(df)

datasets = ['insilico_size10_1','insilico_size10_2','insilico_size10_3','insilico_size10_4','insilico_size10_5']

label_list = []
auroc_list = []

for i in range(1,6):
    for test_statistic in ["auroc","aupr"]:
        main(i, test_statistic)

