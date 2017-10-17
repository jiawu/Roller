import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.append('../pipelines')
import pdb
import Pipelines as pl
import numpy as np
import pandas as pd

td_edge_lists_lag = []

run_params = {}
inf_list = ["Dionesus", "RandomForest", "Lasso"]
inf_list2 = ["Dionesus-td10", "RandomForest-td10", "Lasso-td10", "Dionesus-td15", "RandomForest-td15", "Lasso-td15"]
organisms = ['Ecoli', 'Yeast']

network_results = []
for organism in organisms:
    for network_index in range(1,21):
        method_results = []
        for inf_method in inf_list:
            td_edge_lists_lag.append([])

            run_params['file_path'] = "/projects/p20519/Roller/data/gnw_insilico/network_data/{}/{}-{}_timeseries.tsv".format(organism, organism, network_index)
            run_params['data_folder'] = "/projects/p20519/roller_output/optimizing_window_size/" + str(inf_method) +"/insilizo_size_10_" + str(network_index)
            run_params['td_window'] = 15
            run_params['min_lag'] = 1
            run_params['max_lag'] = 3
            run_params['n_trees'] = 500
            run_params['bootstrap_n'] = 50
            
            roc, pr, tdr = pl.get_td_stats(**run_params)
            result_table = tdr.make_sort_df(tdr.edge_dict, sort_by = 'rank')
            result_table['rank_importance'] = np.arange(len(result_table))
            method_results.append(result_table)

        for i in range(0,3):
          method_results[i] = method_results[i].rename(columns = {"rank_importance":"rank_importance_" + inf_list[i]})
        merged_results = method_results[0].merge(method_results[1], on='regulator-target').merge(method_results[2], on='regulator-target')
        network_results.append(merged_results)

import pickle
pickle.dump(network_results, open("merged_results_td_window_15_all.p","wb"))
#pdb.set_trace()

network_results_td = pickle.load(open("merged_results_td_all.p","rb"))
network_results_td_15 = pickle.load(open("merged_results_td_window_15_all.p","rb"))
network_results = pickle.load(open("merged_results.p_all","rb"))
merged_results = network_results[0]
merged_results_td = network_results_td[0]
merged_results_td = merged_results_td.rename(columns={"rank_importance_Dionesus":"rank_importance_Dionesus_td","rank_importance_Lasso":"rank_importance_Lasso_td","rank_importance_RandomForest":"rank_importance_RandomForest_td"})

merged_results_td_15 = network_results_td_15[0]
merged_results_td_15 = merged_results_td_15.rename(columns={"rank_importance_Dionesus":"rank_importance_Dionesus_td_15","rank_importance_Lasso":"rank_importance_Lasso_td_15","rank_importance_RandomForest":"rank_importance_RandomForest_td_15"})

pdb.set_trace()

mm = merged_results.merge(merged_results_td, on="regulator-target").merge(merged_results_td_15, on="regulator-target")

from sklearn.decomposition import PCA

X = mm.drop('regulator-target',1).T
pca = PCA(n_components = 3)
X_r = pca.fit(X).transform(X)
f = plt.figure(figsize = (10,10) )
axes = f.gca()
inf_list.extend(inf_list2)
for c, i, target_name in zip("rgbkkkkkk", [0, 1, 2,3,4,5, 6,7,8,], inf_list):
      axes.scatter(X_r[i][1], X_r[i][2], c=c, label=target_name)
      axes.annotate(inf_list[i], (X_r[i][1], X_r[i][2]))

f.savefig('pca_plot.ps', format = "ps", bbox_inches='tight')


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

