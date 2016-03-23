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

from lag_identification import get_experiment_list, xcorr_experiments, calc_edge_lag

def get_network_changes(gold_file, pickle_file, edge_str='regulator-target',
                        base_str='rank_importance_RF-td_21', shortener_str='rank_importance_'):
    dg = nx.DiGraph()
    evaluator = Evaluator(gold_file, '\t')
    true_edges = evaluator.gs_flat.tolist()
    dg.add_edges_from(true_edges)
    results_df = pd.read_pickle(pickle_file)
    edges = results_df[edge_str].values
    baseline = results_df[base_str].values

    df = pd.read_csv(data_file, sep="\t")
    gene_list = df.columns.values[1:].tolist()
    experiment_list = get_experiment_list(data_file, 21, 10)
    xcorr_array = xcorr_experiments(experiment_list)
    edge_lags = calc_edge_lag(xcorr_array, gene_list, 0.1, 0.3, timestep=1)
    true_lags = edge_lags[edge_lags['Edge'].isin(true_edges)]

    new_df = pd.DataFrame()
    new_df[edge_str] = edges
    new_df['base_rank'] = baseline
    for column in results_df.columns:
        if column != edge_str and column != base_str:
            short_name = column.replace(shortener_str, "")
            new_df[short_name] = baseline - results_df[column].values

    return new_df

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


if __name__ == "__main__":
    net = 1
    gold_file = "../data/gnw_insilico/network_data/Ecoli/Ecoli-%i_goldstandard.tsv" % net
    data_file = "../data/gnw_insilico/network_data/Ecoli/Ecoli-%i_timeseries.tsv" % net
    pickle_file = "Ecoli_net%i_promotion.pkl" % net
    print(get_network_changes(gold_file, pickle_file).head())

    # Wilcoxon test for each network to see if it changed? - may not be appropriate
    # Hypergeometric test for promoted edges enriched for lagged edges?
        # Calculate for each experimental condition
        # Get true edge lag





"""
Old plotting code

demoted = new_df[new_df['diff'] <= 0]
promoted = new_df[new_df['diff'] > 0]
p_edges = set(promoted['regulator-target'].values)
d_edges = set(demoted['regulator-target'].values)
print(len(p_edges.intersection(set(true_edges))), len(d_edges.intersection(set(true_edges))))
# print(len(d_edges.intersection(set(true_edges))))
# print(len(true_edges))
for row in demoted.iterrows():
    print(row[1])
    sys.exit()
    c = 'r'
    plt.plot([0, 1], [row[1]['base_rank'], row[1][column]], color=c)
for row in promoted.iterrows():
    c = 'b'
    plt.plot([0, 1], [row[1]['base_rank'], row[1][column]], color=c)
plt.xticks([0, 1], ['base_rank', column.replace('rank_importance_', "")])
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()
"""






