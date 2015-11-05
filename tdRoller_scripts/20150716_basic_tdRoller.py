__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'
import pdb
from Swing import Swing
from Swing.util.Evaluator import Evaluator
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

insilico_n = 4
window_width = 18
min_lag = 4
max_lag = 5
n_trees = 10
n_permutes = 10
mse_adjust = False
combine_method = 'mean_mean'
sort_by = 'adj'
file_path = "../data/dream4/insilico_size10_%i_timeseries.tsv"%insilico_n
gene_start_column = 1
time_label = "Time"
separator = "\t"
gene_end = None
current_gold_standard = file_path.replace("timeseries.tsv","goldstandard.tsv")
evaluator = Evaluator(current_gold_standard, '\t')
true_edges = evaluator.gs_flat.tolist()
pd.options.display.float_format = '{:,.5f}'.format

np.random.seed(8)

tdr = Swing(file_path, gene_start_column, gene_end, time_label, separator, min_lag=min_lag, max_lag=max_lag)
#tdr.zscore_all_data()
tdr.set_window(window_width)
tdr.create_windows()
tdr.augment_windows(min_lag=min_lag, max_lag=max_lag)
tdr.fit_windows(n_trees=n_trees, show_progress=False, calc_mse=mse_adjust)
tdr.rank_edges(permutation_n=n_permutes, calc_mse=mse_adjust)
tdr.compile_roller_edges(self_edges=True, calc_mse=mse_adjust)
tdr.make_static_edge_dict(true_edges, lag_method=combine_method)
df2 = tdr.make_sort_df(tdr.edge_dict, sort_by)
df2['Rank'] = np.arange(len(df2))

#print(df2)
roc_dict, pr_dict = tdr.score(df2)
#tdr.plot_scoring(roc_dict, pr_dict)

print df2

gs_ranks = [df2['Rank'][df2['regulator-target'] == edge].values[0] for edge in true_edges]
print zip(true_edges, gs_ranks)
# box_data = [tdr.edge_dict[(edge)]['dataframe']['adj_imp'].values for edge in df2['regulator-target'].values]
print tdr.edge_dict[('G10', 'G7')]['dataframe']

print 'AUROC: ', roc_dict['auroc'][-1]
print 'AUPR: ', pr_dict['aupr'][-1]#+(1-pr_dict['recall'][-1])
print 'F1: ', (2*roc_dict['auroc'][-1]*pr_dict['aupr'][-1]/(roc_dict['auroc'][-1]+pr_dict['aupr'][-1]))

print "Network: ", insilico_n
print 'Window width: ', window_width
print 'Min lag: ', min_lag
print 'Max lag: ', max_lag
print 'trees: ', n_trees
print 'perms: ', n_permutes
print 'mse adjusted: ', mse_adjust
print 'lumping: ', combine_method
print 'sorting by: ', sort_by

# plt.boxplot(box_data)
# plt.show()
