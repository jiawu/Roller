import matplotlib
matplotlib.use('Agg')
from Roller.util.LinePlot import LinePlot
from Roller.util.Analyzer import Analyzer

import pdb

from Roller.tdRoller import tdRoller
from Roller.util.Evaluator import Evaluator
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


def get_td_stats(file_path): 
    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None

    df = pd.read_csv(file_path,sep=separator)
    current_gold_standard = file_path.replace("timeseries.tsv","goldstandard.tsv")
    node_list = df.columns.tolist()
    node_list.pop(0)

    evaluator = Evaluator(current_gold_standard, '\t', node_list=node_list)
    true_edges = evaluator.gs_flat.tolist()
    pd.options.display.float_format = '{:,.5f}'.format

    tdr = tdRoller(file_path, gene_start_column, gene_end, time_label, separator)
    tdr.zscore_all_data()
    tdr.set_window(8)
    tdr.create_windows()
    tdr.augment_windows(min_lag=1, max_lag=4)
    tdr.fit_windows(n_trees=10, show_progress=False)
    tdr.rank_edges(permutation_n=10)
    tdr.compile_roller_edges(self_edges=True)

    tdr.full_edge_list.loc[tdr.full_edge_list.p_value>=0.05, 'Importance'] = 0
    tdr.make_static_edge_dict(true_edges, lag_method='median_median')
    df2 = tdr.make_sort_df(tdr.edge_dict, 'lag')
    print len(df2)
    roc_dict, pr_dict = tdr.score(df2)
    print roc_dict['auroc'][-1]
    print pr_dict['aupr'][-1]#+(1-pr_dict['recall'][-1])
    #tdr.plot_scoring(roc_dict, pr_dict)
    #return((roc_dict['auroc'][-1],pr_dict['aupr'][-1]))
    return((roc_dict, pr_dict))

dataset_list = ['1','2','3','4','5']
for data_index in dataset_list:
    data_folder = "/projects/p20519/roller_output/optimizing_window_size/RandomForest/insilico_size10_" + data_index + "/"

    output_path = "/home/jjw036/insilico_"+data_index+"_"

    target_dataset = "/projects/p20519/Roller/data/dream4/insilico_size10_"+data_index+"_timeseries.tsv"
    roc,pr = get_td_stats(target_dataset)
    #target_dataset = "/projects/p20519/Roller/data/invitro/whitfield_shojaie_timeseries.tsv"

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

    max_results = analyzer.max_width_results
    crag_results = analyzer.min_width_results
    tdr_label = "tdSWING, aupr = " + format(pr['aupr'][-1], '.2f')
    crag_label = "SWING, aupr = " + format(crag_results['aupr'].tolist()[-1], '.2f')
    max_label = "RandomForest, aupr = " + format(max_results['aupr'].tolist()[-1], '.2f')
    lp.plot_window_series(pr['precision'],1, x_values=pr['recall'], label=tdr_label)
    lp.plot_window_series(crag_results['precision'],2, x_values=crag_results['recall'], label=crag_label)
    lp.plot_window_series(max_results['precision'],3, x_values=max_results['recall'], label=max_label)

    ## plot vertical line for the maximum window size
    #lp.plot_vertical_line(analyzer.current_roller.overall_width, color_index, window_size)

    ## print best cragging score
    cragged_window = analyzer.predict_best_window()
    #max_window = analyzer.get_max_window()
    #lp.plot_horizontal_line(cragged_window['auroc'].values, 1, 'best crag')
    lp.add_formatting_auroc()        

    lp.save_plot(output_path, 'aupr')

    #grouped.get_group((2,2)).mean()['aupr']


