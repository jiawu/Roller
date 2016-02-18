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

"""
Script that loads data from a dataframe and generates boxplots

"""

def read_tdr_results(folder_list):
    agg_df = pd.DataFrame()
    for input_folder in folder_list:
        for file_path in os.listdir(input_folder):    
            df = pd.read_csv(input_folder+file_path,sep=',|\t', engine='python')
            agg_df = agg_df.append(df)
    return(agg_df)

def parse_tdr_results(agg_df,test_statistic, datasets):
    label_list = []
    auroc_list = []

    ## Analyze:
      # nonuniform
      # uniform
      # for all networks 1 2 3 4 5
      # parsing for windows = 7, windows = 4

    for dataset in datasets:
        current_df = agg_df[agg_df['file_path'].str.contains(dataset)]
        #current_df = current_df[current_df['iterating_param']=='filter_noisy']


        RF = current_df[(current_df['td_window'] == 21)]
        granger_RF = current_df[(current_df['td_window'] == 20) & (current_df['min_lag']==1) ]
        SWING_RF = current_df[(current_df['td_window'] == 10) & (current_df['min_lag']==1) & (current_df['filter_noisy'] == 0)]
        SWING_RF_filtered = current_df[(current_df['td_window'] == 10) & (current_df['min_lag']==1) & (current_df['filter_noisy'] == 1)]

        comparisons = [RF, granger_RF, SWING_RF, SWING_RF_filtered]
        
        for category in comparisons:
            auroc_list.append(category[test_statistic][0:n_trials].tolist())
        
        label_list.append("RF")
        label_list.append("Granger RF")
        label_list.append("SWING RF")
        label_list.append("SWING RF-Filtered")

    
    return((label_list, auroc_list))

output_path = "/home/jjw036/"

input_folder_list = ["/projects/p20519/roller_output/gnw/RandomForest/"]  
test_statistic = ['aupr', 'auroc']
save_tag = "RF_comparisons_detailed"
n_trials = 100

#datasets = ["_"]
datasets = ["Ecoli-"+str(index)+"_" for index in range(1,6)]
#datasets = ['insilico_size10_1','insilico_size10_2','insilico_size10_3','insilico_size10_4','insilico_size10_5']
agg_df = read_tdr_results(input_folder_list)
with PdfPages(output_path+save_tag+'.pdf') as pdf:
    for test in test_statistic:
        label_list, auroc_list = parse_tdr_results(agg_df,test, datasets)
        bp_data = auroc_list
        bp = BoxPlot()
        bp.plot_box(bp_data, label_list)
        title = save_tag
        bp.add_formatting(title, y_label=test.upper())
        labels = ["Net "+str(x) for x in range(1,6)]
        bp.add_sections(4, labels, offset=0.1)
        
        scoring_scheme = [(x, x+1) for x in range(1,20, 4)]
        scoring_scheme = scoring_scheme + [(x,x+1) for x in range(2,20,4)]
        #scoring_scheme = [(x,x+1) for x in range(0, 19, 2)]
        #print(scoring_scheme)
        tests = bp.sigtest(bp_data, score=scoring_scheme)
        print(tests)
        print(len(bp_data[0]))
        bp.add_significance(tests)
        pdf.savefig(bp.f)
   



#auroc_1 = df['auroc'].values
#auroc_2 = df['auroc'].values

#bp_data = [auroc_1,auroc_2]

#bp = BoxPlot()

#bp.plot_box(bp_data, ['n_trees = 10', 'n_trees = 20'])


#bp.save_plot(output_path, save_tag)

    

    #grouped.get_group((2,2)).mean()['aupr']


