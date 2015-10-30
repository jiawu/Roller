__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

from Swing.tdSwing import tdSwing
from Swing.util.Evaluator import Evaluator
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import pdb

file_path = "../data/invitro/janes_timeseries.tsv"
gene_start_column = 1
time_label = "Time"
separator = "\t"
gene_end = None

df = pd.read_csv(file_path,sep=separator)
current_gold_standard = file_path.replace("timeseries.tsv","goldstandard.tsv")
node_list = df.columns.tolist()
node_list.pop(0)

evaluator = Evaluator(current_gold_standard, '\t', node_list=node_list)
true_edges = evaluator.gs_flat.tolist()
pd.options.display.float_format = '{:,.5f}'.format

np.random.seed(8)

tdr = tdSwing(file_path, gene_start_column, gene_end, time_label, separator)
tdr.zscore_all_data()
tdr.set_window(8)
tdr.create_windows()
tdr.augment_windows(min_lag=1, max_lag=4)
tdr.fit_windows(n_trees=10, show_progress=False)
tdr.rank_edges(permutation_n=10)
tdr.compile_roller_edges(self_edges=True)
pdb.set_trace()

tdr.full_edge_list.loc[tdr.full_edge_list.p_value>=0.05, 'Importance'] = 0
tdr.make_static_edge_dict(true_edges, lag_method='median_median')
df2 = tdr.make_sort_df(tdr.edge_dict, 'lag')
print len(df2)
roc_dict, pr_dict = tdr.score(df2)
print roc_dict['auroc'][-1]
print pr_dict['aupr'][-1]#+(1-pr_dict['recall'][-1])
#tdr.plot_scoring(roc_dict, pr_dict)
