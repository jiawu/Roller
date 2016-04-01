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

    ## Analyze:
      # nonuniform
      # uniform
      # for all networks 1 2 3 4 5
      # parsing for windows = 7, windows = 4
    cum_aurocs = []
    for dataset in datasets:
        auroc_list = []
        current_df = agg_df[agg_df['file_path'].str.contains(dataset)]

        RF = current_df[(current_df['td_window'] == 21) & (current_df['data_folder'].str.contains('RandomForest')) & (current_df['data_folder'].str.contains('yeast'))]
        SWING_RF = current_df[(current_df['td_window'] == 5) & (current_df['min_lag']==0) & (current_df['max_lag'] == 1) & (current_df['data_folder'].str.contains('RandomForest')) & (current_df['data_folder'].str.contains('yeast'))]
        SWING_Lasso = current_df[(current_df['td_window'] == 10) & (current_df['min_lag']==0) & (current_df['max_lag'] == 1) & (current_df['data_folder'].str.contains('RandomForest'))& (current_df['data_folder'].str.contains('yeast') )]
        SWING_Dionesus = current_df[(current_df['td_window'] == 15) & (current_df['min_lag']==0) & (current_df['max_lag'] == 1) & (current_df['data_folder'].str.contains('RandomForest'))& (current_df['data_folder'].str.contains('yeast') )]
        SWING_Community = current_df[(current_df['td_window'] == 20) & (current_df['min_lag']==0) & (current_df['max_lag'] == 1) & (current_df['data_folder'].str.contains('RandomForest'))& (current_df['data_folder'].str.contains('yeast') )]
        #SWING_Community = current_df[(current_df['td_window'] == 18) & (current_df['max_lag'] == 3) & (current_df['data_folder'].str.contains('RandomForest'))]


        RF2 = current_df[(current_df['td_window'] == 21) & (current_df['data_folder'].str.contains('RandomForest'))& (current_df['data_folder'].str.contains('ecoli') )]
        SWING_RF2 = current_df[(current_df['td_window'] == 5) & (current_df['min_lag']==0) & (current_df['max_lag'] == 1) & (current_df['data_folder'].str.contains('RandomForest'))& (current_df['data_folder'].str.contains('ecoli') )]
        SWING_Lasso2 = current_df[(current_df['td_window'] == 10) & (current_df['min_lag']==0) & (current_df['max_lag'] == 1) & (current_df['data_folder'].str.contains('RandomForest'))& (current_df['data_folder'].str.contains('ecoli') )]
        SWING_Dionesus2 = current_df[(current_df['td_window'] == 15) & (current_df['min_lag']==0) & (current_df['max_lag'] == 1) & (current_df['data_folder'].str.contains('RandomForest'))& (current_df['data_folder'].str.contains('ecoli') )]
        SWING_Community2 = current_df[(current_df['td_window'] == 20) & (current_df['min_lag']==0) & (current_df['max_lag'] == 1) & (current_df['data_folder'].str.contains('RandomForest'))& (current_df['data_folder'].str.contains('ecoli') )]

        comparisons = [RF, SWING_RF, SWING_Lasso, SWING_Dionesus, SWING_Community,RF2, SWING_RF2, SWING_Lasso2, SWING_Dionesus2, SWING_Community2 ]

        
        for category in comparisons:
            auroc_list.append(category[test_statistic][0:n_trials].tolist())

        #for each dataset, add up the test statistics. Take the mean of the reference box. Subtract each item from the mean
        reference_box = auroc_list[0]
        ref_mean = np.mean(reference_box)

        for box in auroc_list[0:5]:
            cum_aurocs.append((box-ref_mean)/ref_mean*100)
        
        reference_box2 = auroc_list[5]
        ref_mean2 = np.mean(reference_box2)

        
        for box in auroc_list[5:]:
            cum_aurocs.append((box-ref_mean2)/ref_mean*100)


    final_RF = np.hstack(cum_aurocs[0::10]).tolist()        
    final_SWING_RF = np.hstack(cum_aurocs[1::10]).tolist()        
    final_SWING_Lasso = np.hstack(cum_aurocs[2::10]).tolist()        
    final_SWING_Dionesus = np.hstack(cum_aurocs[3::10]).tolist()        
    final_SWING_Community = np.hstack(cum_aurocs[4::10]).tolist()        
    
    final_RF2 = np.hstack(cum_aurocs[5::10]).tolist()        
    final_SWING_RF2 = np.hstack(cum_aurocs[6::10]).tolist()        
    final_SWING_Lasso2 = np.hstack(cum_aurocs[7::10]).tolist()        
    final_SWING_Dionesus2 = np.hstack(cum_aurocs[8::10]).tolist()        
    final_SWING_Community2 = np.hstack(cum_aurocs[9::10]).tolist()

    final_RF = final_RF + final_RF2
    final_SWING_RF = final_SWING_RF + final_SWING_RF2
    final_SWING_Lasso = final_SWING_Lasso + final_SWING_Lasso2
    final_SWING_Dionesus = final_SWING_Dionesus + final_SWING_Dionesus2
    final_SWING_Community = final_SWING_Community + final_SWING_Community2

    final_auroc_list = [final_RF, final_SWING_RF, final_SWING_Lasso, final_SWING_Dionesus, final_SWING_Community]
        
    label_list.append("RF")
    label_list.append("SWING RF - w=5")
    label_list.append("SWING RF - w=10")
    label_list.append("SWING RF - w=15")
    label_list.append("SWING RF - w=18")
    
    return((label_list, final_auroc_list))

output_path = "/home/jjw036/"

input_folder_list = ["/projects/p20519/roller_output/community/", "/projects/p20519/roller_output/gnw/RandomForest/", "/projects/p20519/roller_output/gnw/Lasso/", "/projects/p20519/roller_output/gnw/Dionesus/"]  
test_statistic = ['aupr', 'auroc']
save_tag = "gnw_window_comparison"
n_trials = 100

#datasets = ["_"]

datasets = ["-"+str(index)+"_" for index in range(1,21)]
#datasets = ['insilico_size10_1','insilico_size10_2','insilico_size10_3','insilico_size10_4','insilico_size10_5']
agg_df = read_tdr_results(input_folder_list)

with PdfPages(output_path+save_tag+'.pdf') as pdf:
    for test in test_statistic:
        label_list, auroc_list = parse_tdr_results(agg_df,test, datasets)

        #200 box plots -> each out being SWING, etc

        #condense into 10, adding up each dataset
        bp_data = auroc_list
        bp = BoxPlot()

        pdb.set_trace()
        bp.plot_box(bp_data, label_list)

        scoring_scheme = [(4,1), (4,2), (4,3),(4,0)]

        tests = bp.sigtest(bp_data, score=scoring_scheme)
        print(tests)
        title = "All Networks"
        bp.add_formatting(title, y_label="% "+ test.upper())
        #labels = ['Yeast', 'E. Coli']
        #bp.add_sections(5, labels, offset=0.25)

        #tests = bp.add_significance(tests, style = 'cascade')
        
        pdf.savefig(bp.f)
   



#auroc_1 = df['auroc'].values
#auroc_2 = df['auroc'].values

#bp_data = [auroc_1,auroc_2]

#bp = BoxPlot()

#bp.plot_box(bp_data, ['n_trees = 10', 'n_trees = 20'])


#bp.save_plot(output_path, save_tag)

    

    #grouped.get_group((2,2)).mean()['aupr']


