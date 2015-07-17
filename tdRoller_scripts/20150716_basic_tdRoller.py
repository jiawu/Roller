__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from Roller.tdRoller import tdRoller
from Roller.util.Evaluator import Evaluator
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

insilico_n = 2
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

tdr = tdRoller(file_path, gene_start_column, gene_end, time_label, separator)
tdr.zscore_all_data()
tdr.set_window(11)
tdr.create_windows()
tdr.augment_windows(min_lag=0, max_lag=None)
tdr.fit_windows(n_trees=10, show_progress=False)
tdr.rank_edges(permutation_n=100)
tdr.compile_roller_edges(self_edges=True)
tdr.full_edge_list = tdr.full_edge_list[tdr.full_edge_list.p_value<0.05]
#tdr.full_edge_list = tdr.full_edge_list[tdr.full_edge_list.Lag>1]
tdr.make_static_edge_dict(true_edges)
df2 = tdr.make_sort_df(tdr.edge_dict, 'lag')
print len(df2)
roc_dict, pr_dict = tdr.score(df2)
print roc_dict['auroc'][-1]+(1-roc_dict['fpr'][-1])
print pr_dict['aupr'][-1]#+(1-pr_dict['recall'][-1])
tdr.plot_scoring(roc_dict, pr_dict)