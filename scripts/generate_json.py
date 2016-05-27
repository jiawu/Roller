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

def get_clist(clusterid, cluster_table):
    clist = cluster_table[cluster_table['__glayCluster'] == clusterid]['name'].tolist()
    return(clist)

def get_tf_list():
    """returns full list of tfs"""
    tf_fp = "/home/jjw036/Roller/data/invitro/omranian_parsed_tf_list.tsv"
    with open(tf_fp, 'r') as tf_file:
        tf_list = tf_file.read().splitlines()
    return(tf_list)

def generate_json(merged_lag,method, CLUSTER=14):
    """
    Generates a json file with the lag information, module information, and parent information.
    """
    # the key is the parent
    # first organize default dict such that each parent has all edges
    parent_map = defaultdict(list)
    for parent,child in merged_lag['index'].tolist():
        parent_map[parent].append(child)
    # then expand the dict into a list of dicts
    json_list = []
    for key in parent_map.keys():
        json_dict = {}
        clusterid = merged_lag[merged_lag['parent']==key]['parent_cluster'].tolist()[0]
        json_dict['name'] = 'Module%d.%s' % (clusterid, key)
        json_dict['imports'] = []
        for value in parent_map[key]:
            child_id = merged_lag[(merged_lag['child']==value) & (merged_lag['parent']==key)]['child_cluster'].tolist()[0]

            edge_info = {}
            edge_info['t_name'] = 'Module%d.%s' % (child_id, value)

            lag =  merged_lag[merged_lag['index'] == (key,value)][method].tolist()[0]
            if math.isnan(lag):
                lag = 0
            edge_info['lag'] = lag
            
            json_dict['imports'].append(edge_info)
        json_list.append(json_dict)

    ### fill in empty dicts for child nodes that do not become parents

    all_child = merged_lag['child'].tolist()
    child_only = list(set(all_child) - set(parent_map.keys()))

    for child in child_only:
        json_dict = {}
        clusterid = merged_lag[merged_lag['child']==child]['child_cluster'].tolist()[0]
        json_dict['name'] = 'Module%d.%s' % (clusterid, child)
        json_dict['imports'] = []
        json_list.append(json_dict)

    json_list = sorted(json_list, key=lambda k: len(k['imports']), reverse=False)
    with open('lagged_network_c'+str(CLUSTER)+'.json', 'w') as fp:
        json.dump(json_list, fp, sort_keys=True)

    # for every parent, append the edge
    # dict with name
    return(True)

def run_subswing(df, td_window=6, min_lag = 0, max_lag = 0, window_type = 'RandomForest'):
    """
    Pass in subnet_dict
    """
    true_edges = df['index'].tolist()
    sub_dict = get_subnetwork_info(df)
    file_path = "/home/jjw036/Roller/data/invitro/omranian_parsed_timeseries.tsv"
    gene_start_column = 1
    gene_end = None
    time_label = "Time"
    separator = "\t"

    tdr = Swing(file_path, gene_start_column, gene_end, time_label, separator, min_lag =min_lag, max_lag = max_lag, window_type = window_type, sub_dict=sub_dict)
    
    # remember the data is already zscored
    #tdr.zscore_all_data()
    
    tdr.set_window(td_window)
    tdr.create_custom_windows(sub_dict['tfs'])
    tdr.optimize_params()
    tdr.crag = False
    tdr.calc_mse = False
    tdr.fit_windows(n_trees=100, show_progress=False, n_jobs=-1)
    tdr.rank_edges(permutation_n=10, n_bootstraps=10)
    tdr.compile_roller_edges(self_edges=False)
    tdr.make_static_edge_dict(true_edges, self_edges=False, lag_method='mean_mean')
    sub_df = tdr.make_sort_df(tdr.edge_dict, sort_by = 'rank')
    sub_df['Rank'] = np.arange(len(sub_df))

    sub_eval = Evaluator(subnet_dict = sub_dict)
    
    pr = sub_eval.calc_pr(sub_df.sort('Rank'))
    roc = sub_eval.calc_roc(sub_df.sort('Rank'))
    return(pr[2].values[-1], roc[2].values[-1])

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
    
def extract_subnetwork(cluster_id, merged_lag, parsed_info, agg_results):    
    
    # get the ranks of all the edges, and only true edges
    # get the name of all the transcription factors
    sub_genes = merged_lag['parent'].unique().tolist() + merged_lag['child'].unique().tolist()
    sub_genes = set(sub_genes)

    tf_list = get_tf_list()
    sub_tfs = sub_genes.intersection(set(tf_list))

    targets = sub_genes
    regulators = tf_list
    evaluator = Evaluator()
    sub_all_edges = tuple(map(tuple,evaluator.possible_edges(np.array(regulators), np.array(list(targets)))))
   
    sub_all_edges = [ x for x in sub_all_edges if x[0] != x[1] ]
    sub_true_edges = merged_lag['index'].tolist()

    sub_stats = { 'edges': sub_all_edges,
                  'true_edges': sub_true_edges,
                  'tfs': sub_tfs,
                  'genes': list(sub_genes)}
    
    
    
    result_group = agg_results.groupby(['data_folder','min_lag','max_lag', 'td_window'])
    
    # first parse the agg_df such that the only files you need are in a dataframe
        # group them in terms of inference method, windowing, and lag
        # groupby
    gs_file = agg_results['file_path'].iloc[0].replace('_timeseries.tsv', '_goldstandard.tsv')
    om_eval = Evaluator(gs_file)
    
    # pass in the sub_group df and the parsed info.
    # return the aupr and auroc
    
    baseline_group = ('/projects/p20519/roller_output/ranks/RandomForest/omranian_', '0', '0', '6')
    swing_group =('/projects/p20519/roller_output/ranks/RandomForest/omranian_', '1', '1', '5') 
    
    baseline_group = result_group.get_group(baseline_group).convert_objects(convert_numeric=True)
    sub_pr1, sub_roc1 = get_group_stats(baseline_group, parsed_info, sub_stats, sub_all_edges)
    
    swing_group = result_group.get_group(swing_group).convert_objects(convert_numeric=True)
    sub_pr2, sub_roc2 = get_group_stats(swing_group, parsed_info, sub_stats, sub_all_edges)
    sub_stats['baseline_pr'] = sub_pr1
    sub_stats['baseline_roc'] = sub_roc1
    sub_stats['swing_pr'] = sub_pr2
    sub_stats['swing_roc'] = sub_roc2
    return(sub_stats)

def get_group_stats(sub_group,parsed_info, sub_stats, sub_all_edges):
    #print(summary)
    sub_ranks = [parsed_info.get(k) for k in sub_group['result_path'].tolist() if k in parsed_info.keys()]
    
    if 'rank_importance' in sub_ranks[0].columns:
        sub_average = ut.average_rank(sub_ranks, 'rank_importance')
    else:
        # This is a community dataframe
        for df in sub_ranks:
            df['mean-rank-dr'] = df[['Dionesus-rank','RandomForest-rank']].mean()
        sub_average = ut.average_rank(sub_ranks, 'mean-rank-dr')
        
    sub_average['regulator-target'] = sub_average['regulator-target'].apply(eval)
    #only get regulator-copy pairs in the gold standard

    sub_eval = Evaluator(subnet_dict = sub_stats)
    
    sub_df = sub_average[sub_average['regulator-target'].isin(sub_all_edges)]
    pr = sub_eval.calc_pr(sub_df.sort('mean-rank'))
    roc = sub_eval.calc_roc(sub_df.sort('mean-rank'))
    return(pr[2].values[-1], roc[2].values[-1])


"""
1. Parse gene ontology to get names of modules
2. Identify clusters that are lagged
3. Determine if clusters have higher AUROC with tdRoller than the baseline community or baseline method
4. Check if there's an enrichment, or if cluster is statistically significant

"""
def main(window_type='RandomForest', CLUSTER=14):
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
    if os.path.isfile('lag_df2_parse_biocyc_3.pkl'):
        lag_df = pd.read_pickle('lag_df2_parse_biocyc_3.pkl')
        edge_df = pd.read_pickle('edge_df2_parse_biocyc_3.pkl')
    else:
        ## Get the lags to associate with the network
        (lag_df, edge_df) = flm.get_true_lags('../data/invitro/omranian_parsed_timeseries.tsv',5,30)
        lag_df['lag_median'] = [np.median(x) for x in lag_df['Lag'].tolist()]
        edge_df['lag_median'] = [np.median(x) for x in edge_df['Lag'].tolist()]
        lag_df.to_pickle('lag_df2_parse_biocyc_3.pkl')
        edge_df.to_pickle('edge_df2_parse_biocyc_3.pkl')
    
    lag_df['lag_counts'] = [len(x) if type(x) is list else 0 for x in lag_df['Lag'].tolist()]
    edge_df['lag_counts'] = [len(x) if type(x) is list else 0 for x in edge_df['Lag'].tolist()]

    clusters = pd.read_csv('../data/invitro/regulon_cluster_assignments'+str(CLUSTER)+'.csv',sep=',')

    new_lag= lag_df.reset_index()

    new_lag[['parent','child']] = new_lag['index'].apply(pd.Series)
    merged_lag = pd.merge(new_lag, clusters[['name','__glayCluster']], how='left', left_on=['parent'], right_on=['name'])
    merged_lag = merged_lag.rename(columns = {'__glayCluster':'parent_cluster'})

    merged_lag = pd.merge(merged_lag, clusters[['name','__glayCluster']], how='left', left_on=['child'], right_on=['name'])
    merged_lag = merged_lag.rename(columns = {'__glayCluster':'child_cluster'})
    
    merged_lag.loc[merged_lag['lag_counts']<7, ['lag_median']] = 0
    merged_lag = merged_lag.fillna(0)
    target_clusters = merged_lag
    grouped_by_cluster = target_clusters.groupby('parent_cluster')
    clusters = target_clusters['parent_cluster'].unique().tolist()

    small_clusters = []
    
    for clusterid in clusters:
        current_group = grouped_by_cluster.get_group(clusterid)
        sub_dict = get_subnetwork_info(current_group)
        if (len(current_group)<10) or (len(sub_dict['tfs']) <3):
            small_clusters.append(int(clusterid))
    print(len(merged_lag['parent_cluster'].unique()))
    merged_lag['parent_cluster'] = merged_lag['parent_cluster'].replace(small_clusters, 99)
    merged_lag['child_cluster'] = merged_lag['child_cluster'].replace(small_clusters, 99)
    ## Read the cluster assignments file and change all the small cluster IDs to 99
    ## save as a parsed cluster assignment file
    clusters = pd.read_csv('../data/invitro/regulon_cluster_assignments'+str(CLUSTER)+'.csv',sep=',')
    parsed_clusters = clusters.copy()
    parsed_clusters['__glayCluster'] = clusters['__glayCluster'].replace(small_clusters,99)
    
    parsed_path = '../data/invitro/regulon_cluster_assignments_parsed'+str(CLUSTER)+'.csv'
    parsed_clusters.to_csv(parsed_path, sep=',', index=False)
    
    print(len(merged_lag['parent_cluster'].unique()))
    merged_lag = merged_lag[merged_lag['parent_cluster'] != 99]
    merged_lag = merged_lag[merged_lag['child_cluster'] != 99]

    generate_json(merged_lag, method = 'lag_median', CLUSTER=CLUSTER)
    
    return(merged_lag)

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        window_type = str(sys.argv[1])
        CLUSTER = int(sys.argv[2])

    else:
        window_type = 'RandomForest'
    main(window_type, CLUSTER=CLUSTER)
        

"""
pdb.set_trace()

result_df = grouped_by_cluster.median()
result_df['count'] = grouped_by_cluster.count()['name_x']

valid_clusters = result_df[result_df['count']>9]['child_cluster'].tolist()


omranian_promotion = pd.read_pickle('/projects/p20519/roller_output/pickles/omranian_promotion.pkl')

promotion_lag = pd.merge(merged_lag, omranian_promotion, how='left', left_on=['index'], right_on=['regulator-target'])

promotion_lag['rank_diff_D'] = promotion_lag['rank_importance_Dionesus-td_6']-promotion_lag['rank_importance_Dionesus-td_4']


eco_go = parse_go()

clusters = pd.read_csv('../data/invitro/regulon_cluster_assignments.csv',sep=',')


obodag = GODag("go-basic.obo")

goeaobj = GOEnrichmentStudy(
        eco_go.keys(), # List of mouse protein-coding genes
        eco_go, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_by']) # defult multipletest correction method
# For each cluster, get a list of genes.
# For each cluster, test the list of the genes for gene ontology enrichment

valid_clusters = clusters['__glayCluster'].unique().tolist()
for clusterid in valid_clusters:
    genes_0 = get_clist(clusterid, clusters)
    goea_results_all = goeaobj.run_study(genes_0)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_by < 0.05]
    print(len(goea_results_sig))
    for result in goea_results_sig:
        print(clusterid,result.name, result.ratio_in_study)


# show an enrichment modules with lag

# divide edges into modules

# create file in json format


# for each node, add list of pathways


# for each pathway, get genes involved in pathway

# for each pathway, get a list of interactions within that pathway
# for each pathway, get a list of interactions related to that pathway

# get list of lags, then for each lag 
"""
