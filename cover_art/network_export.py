import pandas as pd
from Swing import Swing
import networkx as nx
import numpy as np
import sys
import matplotlib.pyplot as plt
from collections import Counter


file_path = "/Volumes/Hephaestus/jfinkle/Documents/Northwestern/MoDyLS/Code/Python/Roller/data/gnw_insilico/high_sampling/Yeast100/Yeast100-1_dream4_timeseries.tsv"
gs = "/Volumes/Hephaestus/jfinkle/Documents/Northwestern/MoDyLS/Code/Python/Roller/data/gnw_insilico/high_sampling/Yeast100/Yeast100-1_goldstandard_signed.tsv"
gene_start_column=1
gene_end = None
time_label = 'Time'
separator = '\t'
tdr = Swing(file_path, gene_start_column, gene_end, time_label, separator)
tdr.zscore_all_data()
edges = pd.read_csv(gs, sep='\t', header=None)
edges.columns = ['parent', 'child', 'sign']
dg = nx.from_pandas_dataframe(edges, 'parent', 'child', 'sign', create_using=nx.DiGraph())       # type: nx.DiGraph
grouped = tdr.norm_data.groupby('Time').mean()
gene_peak = pd.DataFrame(grouped.abs().idxmax(), columns=['max_time'])
gene_peak['max_value'] = [grouped.loc[t, idx] for idx, t in gene_peak['max_time'].iteritems()]

gene_peak['out_degree'] = gene_peak.index.map(dg.out_degree)
gene_peak['in_degree'] = gene_peak.index.map(dg.in_degree)
gene_peak.sort_values(['out_degree','in_degree'], ascending=[False, True], inplace=True)
print(gene_peak[gene_peak['in_degree']==0]) 

# These might yield good results
parent = "G39"
experiment = 55
filter = np.linspace(0, 1000, 21)
data = tdr.norm_data.loc[experiment*len(tdr.time_vec):experiment*len(tdr.time_vec)+len(tdr.time_vec)-1].set_index('Time').loc[filter]
nodes = ['G39', "G7", "G8", 'G9']
plt.plot(data.index, data[nodes])
plt.legend(['G39', "G7", "G8", 'G9'])



print(edges[edges.parent==parent].head())
print(gene_peak.loc[edges[edges.parent==parent]['child'].values])
print(Counter(pd.Series(nx.single_source_shortest_path_length(dg, parent))))
# plt.show()