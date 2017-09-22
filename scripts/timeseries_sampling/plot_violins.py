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
"""
Script that loads data from a dataframe and generates boxplots

"""

def read_tdr_results(folder_list, folder_str):
    agg_df = pd.DataFrame()
    for input_folder in folder_list:
        for file_path in os.listdir(input_folder):
            if folder_str in file_path:
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
