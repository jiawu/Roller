__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import os
import sys
import time
import pandas as pd
import numpy as np
from Evaluator import Evaluator
from sklearn.neighbors import DistanceMetric
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.cluster.hierarchy as sch
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

    #Add the roller objects to the dataset key
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
          dataset_dict[dataset_key] = {'roller_list': [], 'results_frame':None, 'canberra_distance': None, 'auroc_difference': None}

        dataset_dict[dataset_key]['roller_list'].append(roller_obj)

    #Compile results for all rollers of a given dataset into a dataframe
    for key, value in dataset_dict.iteritems():
        current_results = make_window_table(value['roller_list'])
        value['results_frame'] = current_results
        value['canberra_distance'] = calc_canberra(current_results)
        value['auroc_difference'] = calc_auroc_diff(current_results)
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
            #print sorted
            importance = unsorted.sort(['importance', 'p_value'], ascending=[False, True])

            #print current_ranking

            tpr, fpr, auroc = evaluator.calc_roc(sorted)
            #print auroc[-1]
            #print filtered.shape
            tpr, fpr, auroc = evaluator.calc_roc(importance)

            #print auroc[-1]
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
    #path = "../../output/Roller_outputs_RF_moretrees/"
    #roller_dict, roller_list = load_roller_pickles(path)
    #pd.to_pickle(roller_dict, "../../output/results_pickles/Roller_outputs_RF.pickle")

    roller_dict = pd.read_pickle("../../output/results_pickles/Roller_outputs_RF.pickle")

    #print len(roller_dict.keys())
    for df in roller_dict.itervalues():
        current_frame = df['results_frame']
        sort_indices= np.argsort(current_frame['AUROC'].values)[::-1]
        sorted_auroc = current_frame['AUROC'].values[sort_indices]
        sorted_auroc_difference = df['auroc_difference'][sort_indices]
        sorted_canberra = df['canberra_distance'][sort_indices]

        # The first row now corresponds to the window with the highest AUROC.
        # Get the indices for ascending canberra distance
        new_sort = np.argsort(sorted_canberra[0])


        rankings = np.array([list(ranking) for ranking in current_frame['Edge_ranking'].values])
        sort_indices = new_sort
        sorted_rankings = rankings[sort_indices]

        # Compute and plot dendrogram.
        fig = plt.figure()

        #sch.set_link_color_palette(['r','g','m','c', 'k', 'y', 'w', 'brown'])

        axdendro = fig.add_axes([0.06, 0.01, 0.9, 0.05])
        Y = sch.linkage(sorted_rankings, metric='canberra')
        Z = sch.dendrogram(Y, orientation='bottom')
        axdendro.set_xticks([])
        axdendro.set_yticks([])
        axdendro.axis("off")

        # Plot distance matrix.
        #axmatrix = fig.add_axes([dm_left, dm_bottom, dm_width, dm_height])
        axmatrix = fig.add_axes([0.06, 0.07, 0.9, 0.5])
        #index = Z['leaves']
        index = new_sort
        sorted_auroc = current_frame['AUROC'].values[index]

        #D_ordered = rankings[sort_indices] # Reorder
        D_ordered = rankings[index] # Reorder
        #nameIdList_ordered = [nameIdList[ii] for ii in index] #Reorder
        im = axmatrix.matshow(D_ordered.T, aspect='auto')
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])
        im.set_cmap(cm.Blues_r)

        #Colorbar
        axcolor = fig.add_axes([0.02, 0.07, 0.02, 0.5])
        cbar=plt.colorbar(im, cax=axcolor, orientation='vertical')
        axcolor.tick_params(labelsize=12, labeltop=True, labelbottom=True)

        auc_ax = fig.add_axes([0.06, 0.6, 0.9, 0.1])
        auc_ax.bar(range(len(sorted_auroc)), sorted_auroc, width=1)
        auc_ax.set_xlim([0, len(sorted_auroc)])

        width_ax = fig.add_axes([0.06, 0.73, 0.9, 0.1])
        #sorted_w = current_frame['Width'].values[sort_indices]
        sorted_w = current_frame['Width'].values[index]
        width_ax.bar(range(len(sorted_w)), sorted_w, width=1)
        width_ax.set_xlim([0, len(sorted_auroc)])
        width_ax.set_ylim([min(sorted_w), max(sorted_w)])

        start_ax = fig.add_axes([0.06, 0.86, 0.9, 0.1])
        sorted_start = current_frame['Start'].values[index]
        start_ax.bar(range(len(sorted_start)), sorted_start, width=1)
        start_ax.set_xlim([0, len(sorted_start)])
        plt.show()
        sys.exit()

    sys.exit()
    results_frame = make_window_table(roller_list)
    auroc_difference = calc_auroc_diff(results_frame)
    canberra_distance = calc_canberra(results_frame)
    sort_indices= np.argsort(results_frame['AUROC'].values)[::-1]
    sorted_auroc = results_frame['AUROC'].values[sort_indices]
    sorted_auroc_difference = auroc_difference[sort_indices]
    sorted_canberra = canberra_distance[sort_indices]
    for ii, auroc in enumerate(sorted_auroc):
        print auroc
        plt.scatter(sorted_canberra[ii], sorted_auroc_difference[ii])
        plt.show()