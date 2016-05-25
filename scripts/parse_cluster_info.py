from datetime import datetime
import sys
import pandas as pd
import pdb
import find_lagged_modules as flm
import numpy as np
import pickle
from collections import defaultdict
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
import json
import math
from Swing import Swing
from Swing.util import utility_module as ut
sys.path.append('/home/jjw036/Roller/pipelines')
import Pipelines as pl
from Swing.util.Evaluator import Evaluator
import os.path

def parse_go():
    go = pd.read_csv('../data/invitro/gene_ontology.tsv', sep='\t')
    genes_go = go.iloc[:,[2,4]]
    genes_go.columns = ['name','GO_ID']
    genes = genes_go['name'].str.lower().tolist()
    go_id = genes_go['GO_ID'].tolist()
    go_tuple = list(zip(genes,go_id))
    eco_go = defaultdict(list)
    for genes,go_id in go_tuple:
        eco_go[genes].append(go_id)
    
    return(eco_go)
def get_subnetwork_info(df):
    sub_genes = df['parent'].unique().tolist() + df['child'].unique().tolist()
    sub_genes = set(sub_genes)

    tf_list = get_tf_list()
    
    sub_tfs = list(sub_genes.intersection(set(tf_list)))
    
    targets = sub_genes
    regulators = sub_tfs
    evaluator = Evaluator()
    sub_all_edges = tuple(map(tuple,evaluator.possible_edges(np.array(regulators), np.array(list(targets)))))
   
    sub_all_edges = [ x for x in sub_all_edges if x[0] != x[1] ]
    sub_true_edges = df['index'].tolist()
    
    sub_stats = { 'edges': sub_all_edges,
                  'true_edges': sub_true_edges,
                  'tfs': sub_tfs,
                  'genes': list(sub_genes)}
    
    
    return(sub_stats)

def get_clist(clusterid, cluster_table):
    clist = cluster_table[cluster_table['__glayCluster'] == clusterid]['name'].tolist()
    return(clist)

def get_tf_list():
    """returns full list of tfs"""
    tf_fp = "/home/jjw036/Roller/data/invitro/omranian_parsed_tf_list.tsv"
    with open(tf_fp, 'r') as tf_file:
        tf_list = tf_file.read().splitlines()
    return(tf_list)

def main(window_type='RandomForest', lag_thresh=4, fill_na=False, median_thresh=2):
    """
    df = pd.read_csv('../data/invitro/ecocyc_database_export_ver1_02.txt',sep='\t')

    # Parse the gene lists for each pathway

    pathway_gene_list = []
    pathway_gene_string = []
    for idx, row in df.iterrows():
        gene_string = row['Genes of pathway']
        parsed = gene_string.replace(' ','').replace('"','').lower().split('//')
        parsed_str = gene_string.replace(' ','').replace('"','').lower().replace('//',' ')
        
        pathway_gene_list.append(parsed)
        pathway_gene_string.append(parsed_str)
        print(parsed)

    df['parsed_genes_list'] = pathway_gene_list
    df['parsed_genes_str'] = pathway_gene_string
    """
    if os.path.isfile('lag_df2_parse_biocyc_2.pkl'):
        lag_df = pd.read_pickle('lag_df2_parse_biocyc_2.pkl')
        edge_df = pd.read_pickle('edge_df2_parse_biocyc_2.pkl')
    else:
        ## Get the lags to associate with the network
        (lag_df, edge_df) = flm.get_true_lags('../data/invitro/omranian_parsed_timeseries.tsv',6,17)
        lag_df['lag_median'] = [np.median(x) for x in lag_df['Lag'].tolist()]
        edge_df['lag_median'] = [np.median(x) for x in edge_df['Lag'].tolist()]
        lag_df.to_pickle('lag_df2_parse_biocyc_2.pkl')
        edge_df.to_pickle('edge_df2_parse_biocyc_2.pkl')
    
    lag_df['lag_counts'] = [len(x) if type(x) is list else 0 for x in lag_df['Lag'].tolist()]
    edge_df['lag_counts'] = [len(x) if type(x) is list else 0 for x in edge_df['Lag'].tolist()]

    clusters = pd.read_csv('../data/invitro/regulon_cluster_assignments7.csv',sep=',')

    new_lag= lag_df.reset_index()

    new_lag[['parent','child']] = new_lag['index'].apply(pd.Series)
    merged_lag = pd.merge(new_lag, clusters[['name','__glayCluster']], how='left', left_on=['parent'], right_on=['name'])
    merged_lag = merged_lag.rename(columns = {'__glayCluster':'parent_cluster'})

    merged_lag = pd.merge(merged_lag, clusters[['name','__glayCluster']], how='left', left_on=['child'], right_on=['name'])
    merged_lag = merged_lag.rename(columns = {'__glayCluster':'child_cluster'})

    within_clusters = merged_lag[merged_lag['parent_cluster'] == merged_lag['child_cluster']]
    between_clusters = merged_lag[merged_lag['parent_cluster'] != merged_lag['child_cluster']]

    target_clusters = within_clusters
    target_clusters.loc[target_clusters['lag_counts']<lag_thresh, ['lag_median']] = np.nan
    target_clusters.loc[target_clusters['lag_counts']<lag_thresh, ['lag_mean']] = np.nan
    if fill_na:
        target_clusters = target_clusters.fillna(0)
    grouped_by_cluster = target_clusters.groupby('parent_cluster')
    clusters = target_clusters['parent_cluster'].unique().tolist()

    cluster_summary = pd.DataFrame()
    clusters = sorted(clusters)
    for clusterid in clusters:
        current_group = grouped_by_cluster.get_group(clusterid)
        sub_dict = get_subnetwork_info(current_group)
        if (len(current_group) < 10) or (len(sub_dict['tfs']) < 3):
            continue
        else:
            non_na = current_group[current_group['lag_median'] > median_thresh]
            plag = len(non_na)/len(current_group)
            result = {'cluster_id': clusterid, 'lag_median':current_group['lag_median'].mean(), 'lag_mean':current_group['lag_mean'].mean(), 'percent_lagged': plag, 'total_edges':len(current_group)}
            cluster_summary = cluster_summary.append(result,ignore_index=True)
    cluster_summary.to_csv('lag_info_merged_c7_4.csv', sep='\t', index=False)
    return(cluster_summary)



if __name__ == "__main__":
    main()
