__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import os
import sys
import pandas as pd
import numpy as np
from Evaluator import Evaluator
from sklearn.neighbors import DistanceMetric
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

def load_roller_pickles(pickle_path):
    '''
    Load pickles into an analyzable data structures
    :param pickle_path:
    :return:
    '''
    file_list = next(os.walk(pickle_path))[2]
    nfiles = len(file_list)

    obj_list = []
    counter = 0

    dataset_dict = {}

    for filename in file_list:
      current_file_path = pickle_path + filename
      roller_obj = pd.read_pickle(current_file_path)
      attributes = dir(roller_obj)
      if any("file_path" in attribute for attribute in attributes):
        counter += 1
        #print(str(counter) + " out of " + str(nfiles))
        obj_list.append(roller_obj)

        dataset_key = roller_obj.file_path
        if dataset_key not in dataset_dict:
          dataset_dict[dataset_key] = []

        dataset_dict[dataset_key].append(roller_obj)

    #print(len(obj_list))
    #print(dataset_dict.keys())
    return dataset_dict, obj_list

def make_window_table(roller_list):
    df = pd.DataFrame(columns=['Timeseries', 'Width', 'Roller_idx', 'Start', 'Goldstandard', 'AUROC','Edges', 'Edge_ranking'])
    count = 0
    for roller in roller_list:
        current_timeseries = roller.file_path
        current_gold_standard = current_timeseries.replace("timeseries.tsv","goldstandard.tsv")
        current_gold_standard = '../../'+current_gold_standard
        current_width = roller.window_width
        evaluator = Evaluator(current_gold_standard, '\t')
        for idx, window in enumerate(roller.window_list):
            start_time = min(window.raw_data['Time'].unique())
            unsorted = window.results_table[np.isfinite(window.results_table['p_value'])]
            edge_sorted = unsorted.sort(['regulator-target'])
            edge_list = edge_sorted['regulator-target'].values
            current_ranking = edge_sorted['p_value-rank'].values.astype(int)
            current_ranking = tuple(max(current_ranking)-current_ranking+1)
            sorted = unsorted.sort(['p_value'], ascending=[True])
            tpr, fpr, auroc = evaluator.calc_roc(sorted)
            df.loc[count] = [current_timeseries, current_width, idx, start_time, current_gold_standard, auroc[-1],
                             edge_list, current_ranking]
            count+=1
    return df

def calc_canberra(data_frame):
    rankings = np.array([list(ranking) for ranking in data_frame['Edge_ranking'].values])
    dist = DistanceMetric.get_metric('canberra')
    can_dist = dist.pairwise(rankings)
    return can_dist

def calc_auroc_diff(data_frame):
    aurocs = data_frame['AUROC'].values
    indices = range(len(aurocs))
    col_idx, row_idx = np.meshgrid(indices, indices)
    auroc_mat = np.abs(aurocs[row_idx] - aurocs[col_idx])
    return auroc_mat

if __name__ == '__main__':
    path = "../../output/Roller_outputs_RF_moretrees/"
    roller_dict, roller_list = load_roller_pickles(path)
    results_frame = make_window_table(roller_list)

    plt.scatter(results_frame['Width'].values, results_frame['AUROC'].values)
    plt.show()
    sys.exit()

    auroc_difference = calc_auroc_diff(results_frame)
    canberra_distance = calc_canberra(results_frame)
    unique_c = np.tril(canberra_distance).flatten()
    unique_a = np.tril(auroc_difference).flatten()

    mask = (unique_c!=0)*(unique_a!=0)
    non_zero_c = unique_c[mask]
    non_zero_a = unique_a[mask]
    print pearsonr(non_zero_a, non_zero_c)
    plt.scatter(non_zero_c, non_zero_a)
    plt.show()
    sys.exit()

    unique = np.tril(canberra_distance).flatten()
    non_zero = unique[unique!=0]
    print len(non_zero)
    print len(unique)
    plt.hist(non_zero, bins=50)
    plt.show()

    unique = np.tril(auroc_difference).flatten()
    non_zero = unique[unique!=0]
    print len(non_zero)
    print len(unique)
    plt.hist(non_zero, bins=50)
    plt.show()