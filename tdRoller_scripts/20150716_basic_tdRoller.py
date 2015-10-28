__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'
import pdb
from Roller.tdRoller import tdRoller
from Roller.util.Evaluator import Evaluator
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

insilico_n = 3
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
tdr.set_window(15)
tdr.create_windows()
tdr.augment_windows(min_lag=1, max_lag=3)
tdr.fit_windows(n_trees=10, show_progress=False)
tdr.rank_edges(permutation_n=10)
tdr.compile_roller_edges(self_edges=True, mse_adjust=False)
tdr.make_static_edge_dict(true_edges, lag_method='mean_mean')
df2 = tdr.make_sort_df(tdr.edge_dict, 'lag')
df2['Rank'] = np.arange(len(df2))
print df2

gs_ranks = [df2['Rank'][df2['regulator-target'] == edge].values[0] for edge in true_edges]
print zip(true_edges, gs_ranks)
box_data = [tdr.edge_dict[(edge)]['dataframe']['Importance'].values for edge in df2['regulator-target'].values]
print tdr.edge_dict[('G1', 'G8')]['dataframe']
plt.boxplot(box_data)
#plt.xticks(df2['regulator-target'].values)
plt.show()

#print(df2)
roc_dict, pr_dict = tdr.score(df2)

print roc_dict['auroc'][-1]
print pr_dict['aupr'][-1]#+(1-pr_dict['recall'][-1])
print (2*roc_dict['auroc'][-1]*pr_dict['aupr'][-1]/(roc_dict['auroc'][-1]+pr_dict['aupr'][-1]))
#tdr.plot_scoring(roc_dict, pr_dict)
