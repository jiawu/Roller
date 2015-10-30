import matplotlib
matplotlib.use('Agg')
from Swing.util.LinePlot import LinePlot
from Swing.util.Analyzer import Analyzer

import pdb

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import td_wrapper as tdw


data_folder = "/projects/p20519/roller_output/optimizing_window_size/RandomForest/insilico_size10_1/"

output_path = "/home/jjw036/Swing/insilico_size10_1"

target_dataset = "/projects/p20519/Swing/data/dream4/insilico_size10_1_timeseries.tsv"
roc,pr = tdw.get_td_stats(target_dataset)
#target_dataset = "/projects/p20519/Swing/data/invitro/whitfield_shojaie_timeseries.tsv"

#Analyzer computes AUROC/AUPR/Cragging Scores and organizes it in a table

analyzer = Analyzer(data_folder)

#identify the x axis in analyzer
time_vec = analyzer.current_roller.time_vec.tolist()

lp = LinePlot()

lp.set_x_values(time_vec)

my_df = analyzer.overall_df
grouped = my_df.groupby(['window_width','window_index'])

## iterate through window_sizes 
unique_window_sizes = list(set(analyzer.overall_df['window_width'].tolist()))
best_window_values = []
for color_index, window_size in enumerate(unique_window_sizes):
    series_y = []
    
    ## get unique indices
    unique_indices = my_df[my_df['window_width']==window_size]['window_index'].unique()

    unique_indices.sort()
    for index in unique_indices:
        value = grouped.get_group((window_size, index)).mean()['auroc']
        series_y.append(value)
        
        unique_indices = list(unique_indices)
        time_values = [time_vec[x] for x in unique_indices]
    best_window_values.append(max(series_y))

lp.plot_window_series(best_window_values,color_index, window_size,x_values=unique_window_sizes)

## plot vertical line for the maximum window size
lp.plot_vertical_line(analyzer.current_roller.overall_width, color_index, window_size)

## print best cragging score
cragged_window = analyzer.predict_best_window()
lp.plot_horizontal_line(cragged_window['auroc'].values, 1, 'best crag')
lp.add_formatting()        

lp.save_plot(output_path, 'test')

#grouped.get_group((2,2)).mean()['aupr']


