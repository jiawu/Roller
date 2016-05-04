import pandas as pd
import pdb

import Swing.util.lag_identification as lag_id
from Swing.util.Evaluator import Evaluator

def get_true_edges(gold_filename):
    evaluator = Evaluator(gold_filename, '\t')
    edges = evaluator.gs_flat.tolist()
    return edges, evaluator

exp_inter = lag_id.get_experiment_list("../data/invitro/omranian_parsed_timeseries.tsv", timepoints=6, perturbs=9)
pdb.set_trace()
print("Got experiments...")
exp_xcor = lag_id.xcorr_experiments(exp_inter, 1)
print("Calculated x correlation matrices")
gene_list = list(exp_inter[0].columns.values)
print("Calculated edge lag from x correlation matrices")

lags=lag_id.calc_edge_lag(exp_xcor, gene_list, 0.1, 0.8, timestep=1)
lags = lags[lags['Parent'] != lags['Child']]
edge_df = pd.DataFrame(lags['Lag'].values, index=lags['Edge'].values, columns=['Lag'])
print(edge_df)


goldstandard = '../data/invitro/omranian_parsed_goldstandard.tsv'
true_edges, evaluator = get_true_edges(goldstandard)

only_true = edge_df[edge_df.index.isin(true_edges)]
pdb.set_trace()

#combine/parse omranian

# Time Gene Names

#map omranian to strong regulondb gold standard

#use swing on gold standard

#check if edges are lagged in omranian

#convert network into modules

#get lagged modules (percentage of each module that is lagged)

#get functional enrichment analysis of modules


