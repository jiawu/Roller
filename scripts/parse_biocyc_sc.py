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
sys.path.append('../pipelines')
import Pipelines as pl
from Swing.util.Evaluator import Evaluator
import os.path
import Swing.util.lag_identification as lag_id

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
    tf_fp = "/home/jjw036/Roller/data/invitro/marbach_parsed_tf_list.tsv"
    with open(tf_fp, 'r') as tf_file:
        tf_list = tf_file.read().splitlines()
    return(tf_list)

def generate_json(merged_df,method):
    """
    Generates a json file with the lag information, module information, and parent information.
    """
    # the key is the parent
    # first organize default dict such that each parent has all edges
    parent_map = defaultdict(list)
    for parent,child in merged_lag['Edge'].tolist():
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

            lag =  merged_lag[merged_lag['Edge'] == (key,value)][method].tolist()[0]
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
    with open('sc_lagged_network.json', 'w') as fp:
        json.dump(json_list, fp, sort_keys=True)

    # for every parent, append the edge
    # dict with name
    return(True)

def run_subswing(df, td_window=6, min_lag = 0, max_lag = 0, window_type = 'RandomForest'):
    """
    Pass in subnet_dict
    """
    true_edges = df['Edge'].tolist()
    sub_dict = get_subnetwork_info(df)
    file_path = "/home/jjw036/Roller/data/invitro/marbach_parsed_timeseries.tsv"
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
    try:
        sub_all_edges = tuple(map(tuple,evaluator.possible_edges(np.array(regulators), np.array(list(targets)))))
   
        sub_all_edges = [ x for x in sub_all_edges if x[0] != x[1] ]
        sub_true_edges = df['Edge'].tolist()
    
    except IndexError:
        sub_stats= {'tfs': []}
        return(sub_stats)
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
    sub_true_edges = merged_lag['Edge'].tolist()

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
    
    baseline_group = ('/projects/p20519/roller_output/ranks/RandomForest/marbach_', '0', '0', '6')
    swing_group =('/projects/p20519/roller_output/ranks/RandomForest/marbach_', '1', '1', '5') 
    
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

def parse_eco_pathways():
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

def get_dbs():
    pd_databases = ['community_agg_rank.pkl','RandomForest_agg_rank.pkl','Lasso_agg_rank.pkl']
    dict_databases = ['community_rank_data.pkl','RandomForest_rank_data.pkl','Lasso_rank_data.pkl']
    db = zip(pd_databases, dict_databases)
    agg_results = pd.DataFrame()
    parsed_info = {}
    target_fp = 'omranian'
    for pd_d, d_d in db:
        df = pd.read_pickle(pd_d)
        df = df[df['file_path'].str.contains(target_fp)]
        agg_results = agg_results.append(df)
        with open(d_d, 'rb') as infile:
            info = pickle.load(infile)
        
        parsed_info.update( {k: info[k] for k in df['result_path'].tolist() if k in info.keys()})
    return(db)

def main(window_type='RandomForest', CLUSTER=1):
    
    if os.path.isfile('sc_lag_df2_parse_biocyc_4.pkl'):
        lag_df = pd.read_pickle('sc_lag_df2_parse_biocyc_4.pkl')
        edge_df = pd.read_pickle('sc_edge_df2_parse_biocyc_4.pkl')
    else:
        ## Get the lags to associate with the network
        #(lag_df, edge_df) = flm.get_true_lags('../data/invitro/marbach_parsed_timeseries.tsv',6,23, dset='marbach')
        experiments = lag_id.get_experiment_list('../data/invitro/marbach_parsed_timeseries.tsv',6,23)
        signed_edge_list = pd.read_csv('../data/invitro/marbach_signed_parsed_goldstandard.tsv',sep='\t',header=None)
        signed_edge_list.columns=['regulator', 'target', 'signs']
        signed_edge_list['regulator-target'] = tuple(zip(signed_edge_list['regulator'],signed_edge_list['target']))
        genes = list(experiments[0].columns.values)
        lag_df,edge_df = lag_id.calc_edge_lag2(experiments,genes,signed_edge_list, mode='marbach')
        lag_df.to_pickle('sc_lag_df2_parse_biocyc_4.pkl')
        edge_df.to_pickle('sc_edge_df2_parse_biocyc_4.pkl')
        
        #lag_df['lag_median'] = [np.median(x) for x in lag_df['Lag'].tolist()]
        #edge_df['lag_median'] = [np.median(x) for x in edge_df['Lag'].tolist()]
        #lag_df.to_pickle('sc_lag_df2_parse_biocyc_3.pkl')
        #edge_df.to_pickle('sc_edge_df2_parse_biocyc_3.pkl')
    
    #lag_df['lag_counts'] = [len(x) if type(x) is list else 0 for x in lag_df['Lag'].tolist()]
    #edge_df['lag_counts'] = [len(x) if type(x) is list else 0 for x in edge_df['Lag'].tolist()]

    clusters = pd.read_csv('../data/invitro/yeast_cluster_assign'+str(CLUSTER)+'.csv',sep=',')

    new_lag= lag_df.reset_index()

    #new_lag[['parent','child']] = new_lag['index'].apply(pd.Series)
    merged_lag = pd.merge(new_lag, clusters[['name','__glayCluster']], how='left', left_on=['parent'], right_on=['name'])
    merged_lag = merged_lag.rename(columns = {'__glayCluster':'parent_cluster'})

    merged_lag = pd.merge(merged_lag, clusters[['name','__glayCluster']], how='left', left_on=['child'], right_on=['name'])
    merged_lag = merged_lag.rename(columns = {'__glayCluster':'child_cluster'})

    #average_lag_over_network = merged_lag['lag_mean'].mean()
    #std_lag_over_network = merged_lag['lag_mean'].std()

    #zero_lag_edges = merged_lag[merged_lag['lag_mean']<1].count()

    within_clusters = merged_lag[merged_lag['parent_cluster'] == merged_lag['child_cluster']]
    between_clusters = merged_lag[merged_lag['parent_cluster'] != merged_lag['child_cluster']]

    target_clusters = within_clusters

    grouped_by_cluster = target_clusters.groupby('parent_cluster')

    #target_clusters.loc[target_clusters['lag_counts']<3, ['lag_median']] = 0


    clusters = target_clusters['parent_cluster'].unique().tolist()
    cluster_summary = pd.DataFrame()

    clusters.sort()
    for clusterid in clusters:
        print(clusterid, len(clusters))
        current_group = grouped_by_cluster.get_group(clusterid)
        total_edges = len(current_group)
        nan_edges = len(current_group[current_group['Lag'].isnull()])
        lagged_edges = len(current_group[current_group['Lag'] >= 10])
        lagged_edges_2 = len(current_group[(current_group['Lag'] >= 20)])
        sub_dict = get_subnetwork_info(current_group)
        if (len(current_group) < 10) or (len(sub_dict['tfs']) < 3):
            continue
        pr1,roc1 = run_subswing(current_group, td_window = 6, min_lag = 0, max_lag = 0, window_type = window_type)
        pr2,roc2 = run_subswing(current_group, td_window = 5, min_lag = 1, max_lag = 1, window_type = window_type)
        pr3,roc3 = run_subswing(current_group, td_window = 5, min_lag = 0, max_lag = 1, window_type = window_type)
        pr4,roc4 = run_subswing(current_group, td_window = 4, min_lag = 0, max_lag = 2, window_type = window_type)
        pr5,roc5 = run_subswing(current_group, td_window = 4, min_lag = 1, max_lag = 2, window_type = window_type)
        pr6,roc6 = run_subswing(current_group, td_window = 4, min_lag = 2, max_lag = 2, window_type = window_type)
        print('Diff pr:', pr1-pr2)
        print('diff roc:', roc1-roc2)
        print(pr1,roc1,pr2,roc2,pr3,roc3,pr4,roc4,pr5,roc5,pr6,roc6)

        
        print('total_edges: %d, nan_edges: %d, lagged_edges: %d, stringently_lagged_edges: %d' % (total_edges, nan_edges, lagged_edges, lagged_edges_2))
        print('percentage nan_edges: %.2f, lagged_edges: %.2f, stringently_lagged_edges: %.2f' % (nan_edges/total_edges, lagged_edges/total_edges, lagged_edges_2/total_edges))

        #subnet_info = extract_subnetwork(clusterid, merged_lag, parsed_info, agg_results)
        #print('AUPR_Swing: ', subnet_info['swing_pr'])
        #print('AUROC_Swing: ', subnet_info['swing_roc'])
        #print('AUPR diff: ', subnet_info['baseline_pr'] - subnet_info['swing_pr'])
        #print('AUROC diff: ', subnet_info['baseline_roc'] - subnet_info['swing_roc'])

        cluster_result = {  'cluster_id': clusterid,
                            'total_edges':total_edges,
                            'nan_edges':nan_edges,
                            'lagged_edges':lagged_edges,
                            'lagged_edges2':lagged_edges_2,
                            'percent_lagged':lagged_edges/total_edges,
                            'percent_lagged2':lagged_edges_2/total_edges,
                            'baseline_auroc':roc1,
                            'baseline_aupr':pr1,
                            'swing_aupr':pr2,
                            'swing_auroc':roc2,
                            'swing_aupr2':pr3,
                            'swing_auroc2':roc3,
                            'swing_aupr3':pr4,
                            'swing_auroc3':roc4,
                            'swing_aupr4':pr5,
                            'swing_auroc4':roc5,
                            'swing_aupr5':pr6,
                            'swing_auroc5':roc6

                            }
        cluster_summary = cluster_summary.append(cluster_result, ignore_index = True)

    current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')    
    cluster_summary.to_csv('sc_cluster_summary_within_c'+str(CLUSTER)+'_' + current_time + '.csv', header=True, index=False, sep='\t')

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        window_type = str(sys.argv[1])
        CLUSTER = int(sys.argv[2])
    else:
        window_type = 'RandomForest'
        CLUSTER = 1
    main(window_type, CLUSTER = CLUSTER) 
