import pandas as pd
from Swing import Swing
from Swing.util.Evaluator import Evaluator
import numpy as np
import Swing.util.utility_module as Rutil
import re
import sys

import pdb


def get_td_stats_test(**kwargs): 
    kwargs.setdefault('td_window',6)
    kwargs.setdefault('min_lag',0)
    kwargs.setdefault('max_lag',4)
    kwargs.setdefault('n_trees',500)
    kwargs.setdefault('permutation_n',10)
    kwargs.setdefault('lag_method','median_median')
    kwargs.setdefault('sort_by', 'rank')
    kwargs.setdefault('calc_mse',False)
    kwargs.setdefault('window_type',None)

    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None
    file_path = kwargs['file_path']

    df = pd.read_csv(file_path,sep=separator)
    current_gold_standard = file_path.replace("timeseries.tsv","goldstandard.tsv")
    node_list = df.columns.tolist()
    node_list.pop(0)

    evaluator = Evaluator(current_gold_standard, '\t', node_list=node_list)
    true_edges = evaluator.gs_flat.tolist()
    pd.options.display.float_format = '{:,.5f}'.format

    td_types = ["Dionesus", "RandomForest", "Lasso"]
    final_edge_list = []
    for type in td_types:
        tdr = tdRoller(file_path, gene_start_column, gene_end, time_label, separator, window_type=type)

        tdr.zscore_all_data()
        tdr.set_window(kwargs['td_window'])
        tdr.create_windows()
        tdr.augment_windows(min_lag=kwargs['min_lag'], max_lag=kwargs['max_lag'])
        tdr.fit_windows(n_trees=kwargs['n_trees'], show_progress=False, calc_mse = kwargs['calc_mse'])
        tdr.rank_edges(permutation_n=kwargs['permutation_n'], calc_mse = kwargs['calc_mse'])
        tdr.compile_roller_edges(self_edges=True, calc_mse = kwargs['calc_mse'])

        tdr.make_static_edge_dict(true_edges, lag_method=kwargs['lag_method'])
        df2 = tdr.make_sort_df(tdr.edge_dict, sort_by = kwargs['sort_by'])
        final_edge_list.append(df2)
        roc_dict, pr_dict = tdr.score(df2)
        #tdr.plot_scoring(roc_dict, pr_dict)
        df2.to_csv("janes_"+type+".csv", sep=",")
    

    print(roc_dict['auroc'][-1])
    print(pr_dict['aupr'][-1])#+(1-pr_dict['recall'][-1])
    #tdr.plot_scoring(roc_dict, pr_dict)
    return(roc_dict['auroc'][-1],pr_dict['aupr'][-1], tdr)
    #lp.plot_horizontal_line(cragged_window[my_statistic].values, 3, 'best crag')

def get_td_stats(**kwargs):
    kwargs.setdefault('td_window',15)
    kwargs.setdefault('min_lag',1)
    kwargs.setdefault('max_lag',3)
    kwargs.setdefault('n_trees',100)
    kwargs.setdefault('permutation_n',10)
    kwargs.setdefault('lag_method','mean_mean')
    kwargs.setdefault('calc_mse',False)
    kwargs.setdefault('bootstrap_n',10)
    kwargs.setdefault('sort_by', 'rank')
    kwargs.setdefault('filter_noisy',False)
    kwargs.setdefault('alpha',None)
    kwargs.setdefault('window_type', None)

    gene_start_column = 1
    time_label = "Time"
    separator = "\t"
    gene_end = None
    file_path = kwargs['file_path']
    
    df = pd.read_csv(file_path,sep=separator)
    if "high_sampling" in file_path:
        
        prefix = file_path.split('/')[-1].split('_')[0]
        pref = file_path.split('/')[-1]
        folder = file_path.rstrip(pref)
        current_gold_standard = "{}{}_goldstandard.tsv".format(folder,prefix)

    else:    
        current_gold_standard = file_path.replace("timeseries.tsv","goldstandard.tsv")
    node_list = df.columns.tolist()
    node_list.pop(0)
    evaluator = Evaluator(current_gold_standard, '\t', node_list=node_list)
    true_edges = evaluator.gs_flat.tolist()
    pd.options.display.float_format = '{:,.5f}'.format

    if kwargs['window_type'] == None:
        inf_method = ['Lasso', 'RandomForest', 'Dionesus']
        matches=re.findall(r"(?=(" + '|'.join(inf_method)+r"))",kwargs['data_folder'])
        kwargs['window_type'] = matches[0]


    tdr = Swing(file_path, gene_start_column, gene_end, time_label, separator, min_lag = kwargs['min_lag'], max_lag = kwargs['max_lag'], window_type = kwargs['window_type'])
    final_edge_list = []

    # turning off zscoring for omranian
    if 'omranian' not in file_path:
        tdr.zscore_all_data()
    tdr.set_window(kwargs['td_window'])
    if 'cantone' in kwargs['file_path']:
        tf_list = ['CBF1','SWI5','ASH1', 'GAL4', 'GAL80']
    elif 'omranian' in kwargs['file_path']:
        with open('../data/invitro/omranian_parsed_tf_list.tsv','r') as f:
            tf_list = f.read().splitlines()
    elif 'dream5' in kwargs['file_path']:
        with open('../data/dream5/insilico_transcription_factors.tsv','r') as f:
            tf_list = f.read().splitlines()
    elif '100-' in kwargs['file_path']:
        tf_list = ['G%s'%x for x in range(1,101)]
    else:
        tf_list = ['G%s'%x for x in range(1,11)]

    # No custom windows for Gardner data
    make_custom = True
    if 'gardner' in kwargs['file_path'] or 'lahav' in kwargs['file_path']:
        make_custom = False

    if make_custom:
        tdr.create_custom_windows(tf_list)
    else:
        tdr.create_windows()

    if kwargs['filter_noisy']:
        tdr.filter_noisy()
    if kwargs['alpha'] is not None:
        tdr.alpha = kwargs['alpha']
    tdr.optimize_params()
    tdr.crag = False
    tdr.calc_mse = kwargs['calc_mse']

    tdr.fit_windows(n_trees=kwargs['n_trees'], show_progress=False, n_jobs=-1)
    tdr.rank_edges(permutation_n=kwargs['permutation_n'], n_bootstraps=kwargs['bootstrap_n'])

    tdr.compile_roller_edges(self_edges=False)

    # For the Lahav data we need to restrict the edge list to only look at P53 targets
    if 'lahav' in kwargs['file_path']:
        subnet_dict = {'true_edges': evaluator.gs_flat.tolist(), 'edges': evaluator.gs_data['regulator-target'].tolist()}
        evaluator = Evaluator(current_gold_standard, '\t', node_list=node_list, subnet_dict=subnet_dict)
        true_edges = evaluator.gs_flat.tolist()

    tdr.make_static_edge_dict(true_edges, self_edges=False, lag_method=kwargs['lag_method'],
                              full_edge_set=set(evaluator.full_list.tolist()))
    df2 = tdr.make_sort_df(tdr.edge_dict, sort_by=kwargs['sort_by'])
    df2['Rank'] = np.arange(len(df2))

    roc_dict, pr_dict = tdr.score(df2, gold_standard_file=current_gold_standard, evaluator=evaluator)

    print(roc_dict['auroc'][-1])
    print(pr_dict['aupr'][-1])#+(1-pr_dict['recall'][-1])
    return roc_dict['auroc'][-1], pr_dict['aupr'][-1], tdr, df2['regulator-target'].tolist()


def main(data_folder, output_path, target_dataset, my_statistic, td_score):
    #Analyzer computes AUROC/AUPR/Cragging Scores and organizes it in a table

    analyzer = Analyzer(data_folder)

    #identify the x axis in analyzer
    time_vec = analyzer.current_roller.time_vec.tolist()

    lp = LinePlot()


    my_df = analyzer.overall_df
    grouped = my_df.groupby(['window_width','window_index'])

    max_series = []
    min_series = []
    mean_series = []
    ## iterate through window_sizes
    unique_window_sizes = list(set(analyzer.overall_df['window_width'].tolist()))
    lp.set_x_values(unique_window_sizes)
    for color_index, window_size in enumerate(unique_window_sizes):
        series_y = []
        ## get unique indices
        unique_indices = my_df[my_df['window_width']==window_size]['window_index'].unique()
        unique_indices.sort()
        for index in unique_indices:
            value = grouped.get_group((window_size, index)).mean()[my_statistic]
            series_y.append(value)

        ## plot horizontal line for the maximum window size
        if window_size == analyzer.current_roller.overall_width:
            lp.plot_horizontal_line(value, 19, "Status Quo")

        else:
            unique_indices = unique_indices.tolist()
            time_values = [time_vec[x] for x in unique_indices]
            max_series.append(max(series_y))
            mean_series.append(mean(series_y))
            min_series.append(min(series_y))

    #print maximum window
    unique_window_sizes.pop()

    lp.plot_window_series(max_series,0, "Max",x_values=unique_window_sizes)
    lp.plot_window_series(mean_series,1, "Mean",x_values=unique_window_sizes)
    lp.plot_window_series(min_series,2, "Min",x_values=unique_window_sizes)


    ## print best cragging score
    cragged_window = analyzer.predict_best_window()

    lp.plot_horizontal_line(cragged_window[my_statistic].values, 3, 'best crag')
    # plot the tdroller score
    lp.plot_horizontal_line(td_score,17, 'tdRoller')
    lp.add_formatting(y_label=my_statistic)


    lp.save_plot(output_path, '_'+my_statistic)
    #grouped.get_group((2,2)).mean()['aupr']

