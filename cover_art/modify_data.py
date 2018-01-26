import pandas as pd
import numpy as np
import networkx as nx

# Load data
node_data = pd.read_csv('base_nodes.txt', sep='\t')
edge_data = pd.read_csv('base_edges.txt', sep='\t')

dg = nx.from_pandas_dataframe(edge_data, 'parent', 'child', create_using=nx.DiGraph())     # type: nx.DiGraph
center = 0

# Set node info
node_data['depth'] = pd.Series(nx.single_source_shortest_path_length(dg, center))
node_data['degree'] = pd.Series(dg.degree())

# Recalibrate the center
rotations = 4
node_data.loc[0, 'degree'] = int(node_data.loc[0, 'degree'] * rotations)


edge_data['parent_depth'] = node_data.loc[edge_data.parent.values, 'depth'].values
edge_data['child_depth'] = node_data.loc[edge_data.child.values, 'depth'].values
node_data.to_csv('node_data.tsv', sep='\t', index=False)
edge_data.to_csv('edge_list.tsv', sep='\t', index=False)
