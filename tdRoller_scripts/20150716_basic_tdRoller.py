__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'
import pdb
from Swing import Swing
from Swing.util.Evaluator import Evaluator
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

pd.set_option('display.width', 1000)
insilico_n = 3
window_width = 20
min_lag = 0
max_lag = 1
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

tdr = Swing(file_path, gene_start_column, gene_end, time_label, separator, min_lag=min_lag, max_lag=max_lag,
            window_type='Lasso')
tdr.zscore_all_data()
tdr.set_window(window_width)
tdr.create_windows()
tdr.optimize_params()
tdr.fit_windows(n_trees=n_trees, show_progress=True, calc_mse=mse_adjust, n_jobs=1, crag=False)
tdr.rank_edges(permutation_n=n_permutes, calc_mse=mse_adjust, n_bootstraps=10)
tdr.compile_roller_edges(self_edges=False, calc_mse=mse_adjust)
tdr.make_static_edge_dict(true_edges, self_edges=False, lag_method=combine_method)
df2 = tdr.make_sort_df(tdr.edge_dict, sort_by)
df2['Rank'] = np.arange(len(df2))

roc_dict, pr_dict = tdr.score(df2)

print df2

gs_ranks = [df2['Rank'][df2['regulator-target'] == edge].values[0] for edge in true_edges]
print zip(true_edges, gs_ranks)

print "Network: ", insilico_n
print 'Window width: ', window_width
print 'Min lag: ', min_lag
print 'Max lag: ', max_lag
print 'trees: ', n_trees
print 'perms: ', n_permutes
print 'mse adjusted: ', mse_adjust
print 'lumping: ', combine_method
print 'sorting by: ', sort_by

print 'AUROC: ', roc_dict['auroc'][-1]
print 'AUPR: ', pr_dict['aupr'][-1]
print 'F1: ', (2*roc_dict['auroc'][-1]*pr_dict['aupr'][-1]/(roc_dict['auroc'][-1]+pr_dict['aupr'][-1]))
