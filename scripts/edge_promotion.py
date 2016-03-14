__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import sys
import pandas as pd
import numpy as np
import networkx as nx
from Swing.util.Evaluator import Evaluator
from nxpd import draw
from nxpd import nxpdParams
import matplotlib.pyplot as plt

#NOTES FROM JIA
"""
Some columns are labeled as RF-td_2, RF-td10, RF-ml_1, etc.

td_X refers to a window size of X. For all these networks, the rest of the parameters are the same:
min_lag is 1, max lag is 3, except RF-td 21, which is just regular random forest.

ml_X refers to a set of min_lags/max_lags for each X.
if X = 0; min_lag = 0, max_lag = 1
elif X = 1; min_lag = 0, max_lag = 2
elif X = 2; min_lag = 0, max_lag = 3
elif X = 3; min_lag = 1, max_lag = 2
elif X = 4; min_lag = 1, max_lag = 4
elif X = 5; min_lag = 2, max_lag = 3
"""

net = 1
gold_file = "../data/gnw_insilico/network_data/Ecoli/Ecoli-%i_goldstandard.tsv" % net
pickle_file = "Ecoli_net%i_promotion.pkl" % net
dg = nx.DiGraph()
evaluator = Evaluator(gold_file, '\t')
true_edges = evaluator.gs_flat.tolist()
dg.add_edges_from(true_edges)
# draw(dg, layout='neato')
a = pd.read_pickle(pickle_file)
edges = a['regulator-target'].values
baseline = a['rank_importance_RF-td_21'].values
for column in a.columns:
    if column != 'regulator-target' and column != 'rank_importance_RF-td_21':
        new_df = pd.DataFrame()
        new_df['regulator-target'] = edges
        new_df['base_rank'] = baseline
        new_df[column] = a[column].values
        new_df['diff'] = baseline - new_df[column]
        demoted = new_df[new_df['diff'] <= 0]
        promoted = new_df[new_df['diff'] > 0]
        for row in demoted.iterrows():
            c = 'r'
            plt.plot([0, 1], [row[1]['base_rank'], row[1][column]], color=c)
        for row in promoted.iterrows():
            print(row[1])
            c = 'b'
            plt.plot([0, 1], [row[1]['base_rank'], row[1][column]], color=c)
        plt.xticks([])
        plt.gca().invert_yaxis()
        plt.show()
        sys.exit()




