import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


import pandas as pd
import pdb
from parse_biocyc import get_subnetwork_info
import itertools
from Swing import Swing

CLUSTER = 26
lag_df = pd.read_pickle('../pickles/lag_df2_parse_biocyc_6.pkl')
edge_df = pd.read_pickle('../pickles/edge_df2_parse_biocyc_6.pkl')

clusters = pd.read_csv('../../data/invitro/regulon_cluster_assignments'+str(CLUSTER)+'.csv',sep=',')
new_lag= lag_df.reset_index()

merged_lag = pd.merge(new_lag, clusters[['name','__glayCluster']], how='left', left_on=['parent'], right_on=['name'])

merged_lag = merged_lag.rename(columns = {'__glayCluster':'parent_cluster'})

merged_lag = pd.merge(merged_lag, clusters[['name','__glayCluster']], how='left',left_on=['child'], right_on=['name'])

merged_lag = merged_lag.rename(columns = {'__glayCluster':'child_cluster'})

within_clusters = merged_lag[merged_lag['parent_cluster'] == merged_lag['child_cluster']]

between_clusters = merged_lag[merged_lag['parent_cluster'] != merged_lag['child_cluster']]


# Get lagged TFs and sort TFs by lag frequency (given that they have more than 10 targets)

lag_freq = merged_lag.groupby('parent')['Lag'].value_counts()
lag_totals = merged_lag['parent'].value_counts()
non_zero_lag = merged_lag[(merged_lag['Lag'] > 0) & (merged_lag['Lag'] < 30)].groupby('parent')['Lag'].count()

TF_lag_stats = pd.DataFrame()
for TF,n_edges in lag_totals.iteritems():
    if TF in non_zero_lag:
      info = {'TF': TF, 'n_edges': n_edges, 'lag_percent': non_zero_lag[TF]/n_edges, 'lag_bins': lag_freq[TF], '# of lagged_edges': non_zero_lag[TF]}
    else:
      info = {'TF': TF, 'n_edges': n_edges, 'lag_percent': 0, 'lag_bins': lag_freq[TF], '# of lagged edges': 0}

    TF_lag_stats = TF_lag_stats.append(info, ignore_index = True)

TF_lag_stats[TF_lag_stats['n_edges'] > 10][['TF', 'lag_percent', 'n_edges', '# of lagged_edges']].to_csv('TF_lag_distribution.csv')


target_clusters = merged_lag
clusters = target_clusters['parent_cluster'].unique().tolist()
valid_clusters = []

grouped_by_cluster = target_clusters.groupby('parent_cluster')

for clusterid in clusters:
    current_group = grouped_by_cluster.get_group(clusterid)
    total_edges = len(current_group)
    sub_dict = get_subnetwork_info(current_group)
    #if (len(current_group) < 10) or (len(sub_dict['tfs']) < 3):
    #    continue
    #else:
    if len(sub_dict['tfs']) < 2:
        continue
    else:
        valid_clusters.append(clusterid)

cluster_lag_stats = pd.DataFrame()
# enumerate all combinations of connections between clusters
for cluster_id in valid_clusters:
    edges = merged_lag[(merged_lag['parent_cluster'] == cluster_id) | (merged_lag['child_cluster'] == cluster_id)]
    n_edges = len(edges)
    lag_bins = edges['Lag'].value_counts()
    n_lags = lag_bins.reset_index(name="count").query("index > 0.0")['count'].values.sum()
    if cluster_id in [43,25]:
        pdb.set_trace()
    if n_lags > 0:
        percent_lags = n_lags/n_edges
    else:
        percent_lags = 0

    info = {'clusterid': cluster_id, 'n_edges': n_edges, 'lag_percent': percent_lags, 'lag_bins': lag_bins, '# of lagged_edges': n_lags}

    cluster_lag_stats = cluster_lag_stats.append(info, ignore_index = True)

cluster_lag_stats[cluster_lag_stats['n_edges'] > 10][['clusterid', 'lag_percent', 'n_edges', '# of lagged_edges']].to_csv('cluster_lag_distribution.csv')
    
"""
for x in itertools.combinations(valid_clusters, 2):
    parent = x[0]
    child = x[1]
    edges= merged_lag[(merged_lag['parent_cluster'] == parent) & (merged_lag['child_cluster'] == child)]

    parent = x[1]
    child = x[0]
    edges = edges.append(merged_lag[(merged_lag['parent_cluster'] == parent) & (merged_lag['child_cluster'] == child)])

    n_edges = len(edges)
    n_lags = edges['Lag'].value_counts()
    percent_lags = n_lags/n_edges

    # get all edges with connections between the two modules
    if (parent == 22) or (child == 22):
        print(edges)
    # get lag distribution
"""

### Print spark lines for lagged transcriptional regulators

# function accepts single target edge. prints each time-course

TF = 'idnr'
lag = 20
interactions = merged_lag[(merged_lag['parent'] == TF) & (merged_lag['Lag'] == lag)]
target_interaction = ('idnr', 'idnd')
parent = 'kdgr'
child = 'edd'


file_path = "/home/jjw036/Roller/data/invitro/iomranian_parsed_timeseries.tsv"
gene_start_column = 1
gene_end = None
time_label = "Time"
separator = "\t"
min_lag = 0
max_lag = 0

window_types = ['Dionesus']
final_edge_list = []

for window_type in window_types:
    tdr = Swing(file_path, gene_start_column , gene_end, time_label, separator, min_lag = min_lag, max_lag = max_lag, window_type = window_type)

    plt.subplot(2,1,1)
    parsed_raw_data = tdr.raw_data
    #time = parsed_raw_data['Time']
    time = [x*10 for x in range(len(parsed_raw_data))]
    plt.plot(time, parsed_raw_data[parent], label = parent)
    plt.plot(time, parsed_raw_data[child], label = child)
    plt.legend(loc = 'best')
    plt.subplot(2,1,2)


    plt.scatter(parsed_raw_data[child], parsed_raw_data[parent], color = 'red')
    plt.scatter(parsed_raw_data[child], parsed_raw_data[parent].shift(-1), color = 'blue')
    plt.savefig('{}.jpg'.format(parent+child))
