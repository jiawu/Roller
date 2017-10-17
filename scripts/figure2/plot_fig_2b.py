import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.append('../../pipelines')
import pdb
import Pipelines as pl
import numpy as np
import pandas as pd
from Swing.util.Evaluator import Evaluator

def parse(method_string):
    
    min_lag = 1
    max_lag = 3
    td_window = 10
    preinf_method = method_string.split('-')[0]
    inf_method = inf_method_key[preinf_method]

    misc = method_string.split('-')[1]

    if "td" in misc:
        td_window = int(misc.split('_')[1])
        if td_window == 21:
            min_lag = 0
            max_lag = 0

    elif "ml" in misc:
        case = int(misc.split('_')[1])
        if case == 0:
            min_lag = 0
            max_lag = 1
        elif case == 1:
            min_lag = 0
            max_lag = 2
        elif case == 2:
            min_lag = 0
            max_lag = 3

        elif case == 3:
            min_lag = 1
            max_lag = 2

        elif case == 4:
            min_lag = 1
            max_lag = 3

        elif case == 5:
            min_lag = 2
            max_lag = 3    

    return(inf_method, td_window, min_lag, max_lag)

td_edge_lists_lag = []

run_params = {}
inf_list = ["Dionesus", "RandomForest", "Lasso"]
inf_list2 = ["Dionesus-td10", "RandomForest-td10", "Lasso-td10", "Dionesus-td15", "RandomForest-td15", "Lasso-td15"]

methods_of_interest = ['RF-td_5', 'RF-td_10', 'RF-td_15','RF-td_21', 'Dionesus-td_5', 'Dionesus-td_10', 'Dionesus-td_15', 'Dionesus-td_21', 'Lasso-td_5', 'Lasso-td_10', 'Lasso-td_15', 'Lasso-td_21', 'RF-ml_0', 'RF-ml_1', 'RF-ml_2', 'RF-ml_3', 'RF-ml_4', 'RF-ml_5', 'Dionesus-ml_1', 'Dionesus-ml_2', 'Dionesus-ml_3', 'Dionesus-ml_4', 'Dionesus-ml_5', 'Lasso-ml_0', 'Lasso-ml_1', 'Lasso-ml_2', 'Lasso-ml_3', 'Lasso-ml_4', 'Lasso-ml_5']

inf_method_key = {'RF':'RandomForest', 'Dionesus': 'Dionesus', 'Lasso':'Lasso'}

organisms = ['Ecoli', 'Yeast']

organism_results = []
for organism in organisms:    
    network_results = []
    for network_index in range(1,21):
        method_results = []
        for method in methods_of_interest:
        
            inf_method, td_window, min_lag, max_lag = parse(method)
            td_edge_lists_lag.append([])

            run_params['file_path'] = "/home/jjw036/Roller/data/gnw_insilico/network_data/"+ organism + "/" + organism + "-" + str(network_index) + "_timeseries.tsv"
            run_params['data_folder'] = "/projects/p20519/roller_output/gnw/" + str(inf_method) +"/" + organism + "-insilico_size10_" + str(network_index)
            run_params['td_window'] = td_window
            run_params['min_lag'] = min_lag
            run_params['max_lag'] = max_lag
            run_params['n_trees'] = 500
            run_params['bootstrap_n'] = 50
            run_params['filter_noisy'] = False
            print(run_params)

            roc, pr, tdr, _ = pl.get_td_stats(**run_params)
            result_table = tdr.make_sort_df(tdr.edge_dict, sort_by = 'rank')
            result_table['rank_importance'] = np.arange(len(result_table))
            method_results.append(result_table)

        for i in range(0,29):
            method_results[i] = method_results[i].rename(columns = {"rank_importance":"rank_importance_" + methods_of_interest[i]})
        current_results = method_results[0]
        for i, frame in enumerate(method_results[1:], 2):
            current_results = current_results.merge(frame, on='regulator-target')
        #merged_results = method_results[0].merge(method_results[1], on='regulator-target').merge(method_results[2], on='regulator-target')
        network_results.append(current_results)
    organism_results.append(network_results)

import pickle
pickle.dump(organism_results, open("all_organism_results.p","wb"))
mm = pickle.load(open("all_organism_results.p","rb"))
mm = mm[0][0]
"""
network_results_td = pickle.load(open("merged_results_td.p","rb"))
network_results_td_15 = pickle.load(open("merged_results_td_window_15.p","rb"))
network_results = pickle.load(open("merged_results.p","rb"))
merged_results = network_results[0]
merged_results_td = network_results_td[0]
merged_results_td = merged_results_td.rename(columns={"rank_importance_Dionesus":"rank_importance_Dionesus_td","rank_importance_Lasso":"rank_importance_Lasso_td","rank_importance_RandomForest":"rank_importance_RandomForest_td"})

merged_results_td_15 = network_results_td_15[0]
merged_results_td_15 = merged_results_td_15.rename(columns={"rank_importance_Dionesus":"rank_importance_Dionesus_td_15","rank_importance_Lasso":"rank_importance_Lasso_td_15","rank_importance_RandomForest":"rank_importance_RandomForest_td_15"})

pdb.set_trace()

mm = merged_results.merge(merged_results_td, on="regulator-target").merge(merged_results_td_15, on="regulator-target")
"""

from sklearn.decomposition import PCA

#mm = current_results
X = mm.drop('regulator-target',1).T
pca = PCA(n_components = 3)
X_r = pca.fit(X).transform(X)
f = plt.figure(figsize = (10,10) )
axes = f.gca()
#inf_list.extend(inf_list2)
#methods_of_interest = ['RF-td_5', 'RF-td_10', 'RF-td_15','RF-td_21', 'Dionesus-td_5', 'Dionesus-td_10', 'Dionesus-td_15', 'Dionesus-td_21', 'Lasso-td_5', 'Lasso-td_10', 'Lasso-td_15', 'Lasso-td_21', 'RF-ml_0', 'RF-ml_1', 'RF-ml_2', 'RF-ml_3', 'RF-ml_4', 'RF-ml_5', 'Dionesus-ml_1', 'Dionesus-ml_2', 'Dionesus-ml_3', 'Dionesus-ml_4', 'Dionesus-ml_5', 'Lasso-ml_0', 'Lasso-ml_1', 'Lasso-ml_2', 'Lasso-ml_3', 'Lasso-ml_4', 'Lasso-ml_5']

indices = [x for x in range(0,len(methods_of_interest))]
for c, i, target_name in zip("bbbryyyrgggryyyrggggggbbbbbbbbb", indices, methods_of_interest):
      axes.scatter(X_r[i][1], X_r[i][2], c=c, label=target_name)
      axes.annotate(methods_of_interest[i], (X_r[i][1], X_r[i][2]))

f.savefig('all_organisms_PCA_plot_pc23.ps', format = "ps", bbox_inches='tight')

for c, i, target_name in zip("bbbryyyrgggryyyrggggggbbbbbbbbb", indices, methods_of_interest):
      axes.scatter(X_r[i][0], X_r[i][1], c=c, label=target_name)
      axes.annotate(methods_of_interest[i], (X_r[i][1], X_r[i][2]))
f.savefig('all_organisms_PCA_plot_pc12.ps', format = "ps", bbox_inches='tight')


# for each inference method, get results
# list of things to compare:
#  DIONESUS
#  LASSO
#  RF
#  td dionesus, lag = 1
#  td lasso
#  td RF
#  lag = 2
#  lag = 3
