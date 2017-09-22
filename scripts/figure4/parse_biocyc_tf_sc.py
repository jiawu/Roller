from datetime import datetime
import sys
import pandas as pd
import pdb
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
import Swing.util.lag_identification as lag_id
import Swing.util.utility_module as Rutil
import networkx as nx
from nxpd import draw


def parse_go():
    go = pd.read_csv('../../data/invitro/gene_ontology.tsv', sep='\t')
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

def run_subswing(df, td_window=6, min_lag = 0, max_lag = 0, window_type = 'RandomForest', clusterid = None, output_fn = None, tfs = None, alpha = None, pcs = None, n_trees = 1000):
    """
    Pass in subnet_dict
    """
    #true_edges = df['Edge'].tolist()
    sub_dict = get_subnetwork_info(df, sub_tfs = tfs)
    true_edges = sub_dict['true_edges']
    #
    #sub_dict['tfs'] = ['gadx','torr','cusr', 'gntr', 'iscr']
    #sub_dict['tfs'] = [clusterid] + sub_dict['tfs'][1:5]
    #sub_dict['tfs'] = list(set(sub_dict['tfs']))
    sub_eval = Evaluator(subnet_dict = sub_dict)
    
    file_path = "/home/jjw036/Roller/data/invitro/imarbach_parsed_timeseries.tsv"
    gene_start_column = 1
    gene_end = None
    time_label = "Time"
    separator = "\t"

    #window_types = ['Dionesus']
    #window_types = ['Lasso']
    window_types = ['RandomForest']
    final_edge_list = []
    for window_type in window_types:
        tdr = Swing(file_path, gene_start_column, gene_end, time_label, separator, min_lag =min_lag, max_lag = max_lag, window_type = window_type, sub_dict=sub_dict)
        
        # remember the data is already zscored
        #tdr.zscore_all_data()
        
        tdr.set_window(td_window)
        tdr.create_custom_windows(sub_dict['tfs'])
        
        if not alpha:
            tdr.optimize_params()
        tdr.crag = False
        tdr.calc_mse = False
        tdr.fit_windows(pcs=pcs, alpha=alpha, n_trees=n_trees, show_progress=False, n_jobs=-1)
        tdr.rank_edges(permutation_n=10, n_bootstraps=10)
        tdr.compile_roller_edges(self_edges=False)
        tdr.make_static_edge_dict(true_edges, self_edges=False, lag_method='mean_mean')
        sub_df = tdr.make_sort_df(tdr.edge_dict, sort_by = 'rank')
        sub_df['Rank'] = np.arange(len(sub_df))
    
        pr = sub_eval.calc_pr(sub_df.sort('Rank'))
        roc = sub_eval.calc_roc(sub_df.sort('Rank'))
        print(window_type,td_window,roc[2].values[-1],pr[2].values[-1])

        final_edge_list.append(sub_df)

    averaged_rank_data = Rutil.average_rank(final_edge_list,'Rank')
    col_names = averaged_rank_data.columns.tolist()
    for i in range(len(window_types)):
        col_names[i] = window_types[i]+'-rank'

    averaged_rank_data.columns = col_names
    averaged_rank_data.sort('mean-rank', inplace=True)
    pr = sub_eval.calc_pr(averaged_rank_data.sort('mean-rank'))
    roc = sub_eval.calc_roc(averaged_rank_data.sort('mean-rank'))
    print('community',td_window,roc[2].values[-1],pr[2].values[-1])
    
    
    #moduleID, source, target, window_size, min_lag, max_lag, true or false, lag or not, lag time, total number of edges in module
    
    sub_df = averaged_rank_data
    
    sub_df['tp'] = sub_df['regulator-target'].isin(sub_eval.gs_flat)
    sub_df=sub_df.merge(df,how = 'outer', right_on='Edge',left_on='regulator-target')
    sub_df['Source'] = [x[0] for x in sub_df['regulator-target'].tolist()]
    sub_df['Target'] = [x[1] for x in sub_df['regulator-target'].tolist()]
    #sub_df = sub_df[(sub_df['parent_cluster'] == sub_df['child_cluster']) | sub_df['parent_cluster'].isnull()]
    sub_df['moduleID'] = clusterid
    sub_df['window_size'] = td_window
    sub_df['min_lag'] = min_lag
    sub_df['max_lag'] = max_lag
    sub_df['total_edges'] = len(df)
    sub_df['pr'] = pr[2].values[-1]
    sub_df['roc'] = roc[2].values[-1]
    sub_df = sub_df.sort('mean-rank')
    
    #sub_df = sub_df.iloc[:len(df)]

    if os.path.isfile(output_fn):
        with open(output_fn,'a') as output:
            sub_df.to_csv(output, header=False, index=False, sep='\t')
    else:
        with open(output_fn,'a') as output:
            sub_df.to_csv(output, header=True, index=False, sep='\t')

    return(pr[2].values[-1], roc[2].values[-1])

def get_subnetwork_info(df, sub_tfs = None):
    sub_genes = df['parent'].unique().tolist() + df['child'].unique().tolist()
    sub_genes = set(sub_genes)

    tf_list = get_tf_list()
    
    if sub_tfs is None:    
        sub_tfs = list(sub_genes.intersection(set(tf_list)))
    #sub_tfs = ['csgd', 'gntr', 'torr']
    #pdb.set_trace()
    
    targets = sub_genes
    regulators = sub_tfs
    evaluator = Evaluator()
    sub_all_edges = tuple(map(tuple,evaluator.possible_edges(np.array(regulators), np.array(list(targets)))))
   
    sub_all_edges = [ x for x in sub_all_edges if x[0] != x[1] ]
    sub_true_edges = df['Edge'].tolist()
    
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


"""
1. Parse gene ontology to get names of modules
2. Identify clusters that are lagged
3. Determine if clusters have higher AUROC with tdRoller than the baseline community or baseline method
4. Check if there's an enrichment, or if cluster is statistically significant

"""
def main(window_type='RandomForest', CLUSTER=4):
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
    if os.path.isfile('../pickles/sc_lag_df2_parse_biocyc_4.pkl'):
        lag_df = pd.read_pickle('../pickles/sc_lag_df2_parse_biocyc_4.pkl')
        edge_df = pd.read_pickle('../pickles/sc_edge_df2_parse_biocyc_4.pkl')
    else:
        experiments = lag_id.get_experiment_list('../../data/invitro/marbach_parsed_timeseries.tsv',6,23)
        signed_edge_list = pd.read_csv('../../data/invitro/marbach_signed_parsed_goldstandard.tsv',sep='\t',header=None)
        signed_edge_list.columns=['regulator', 'target', 'signs']
        signed_edge_list['regulator-target'] = tuple(zip(signed_edge_list['regulator'],signed_edge_list['target']))
        genes = list(experiments[0].columns.values)
        lag_df,edge_df = lag_id.calc_edge_lag2(experiments,genes,signed_edge_list, mode='marbach')
        ## Get the lags to associate with the network
        #(lag_df, edge_df) = flm.get_true_lags('../../data/invitro/omranian_parsed_timeseries.tsv',5,26)
        #lag_df['lag_median'] = [np.median(x) for x in lag_df['Lag'].tolist()]
        #edge_df['lag_median'] = [np.median(x) for x in edge_df['Lag'].tolist()]
        lag_df.to_pickle('../pickles/sc_lag_df2_parse_biocyc_4.pkl')
        edge_df.to_pickle('../pickles/sc_edge_df2_parse_biocyc_4.pkl')
    
    #lag_df['lag_counts'] = [len(x) if type(x) is list else 0 for x in lag_df['Lag'].tolist()]
    #edge_df['lag_counts'] = [len(x) if type(x) is list else 0 for x in edge_df['Lag'].tolist()]

    clusters = pd.read_csv('../../data/invitro/yeast_cluster_assign'+str(CLUSTER)+'.csv',sep=',')

    new_lag= lag_df.reset_index()
    #new_lag[['parent','child']] = new_lag['index'].apply(pd.Series)
    merged_lag = pd.merge(new_lag, clusters[['name','__glayCluster']], how='left', left_on=['parent'], right_on=['name'])
    merged_lag = merged_lag.rename(columns = {'__glayCluster':'parent_cluster'})

    merged_lag = pd.merge(merged_lag, clusters[['name','__glayCluster']], how='left', left_on=['child'], right_on=['name'])
    merged_lag = merged_lag.rename(columns = {'__glayCluster':'child_cluster'})

    #generate_json(merged_lag, method = 'lag_median')
    
    #average_lag_over_network = merged_lag['lag_mean'].mean()
    #std_lag_over_network = merged_lag['lag_mean'].std()

    #zero_lag_edges = merged_lag[merged_lag['lag_mean']<1].count()

    within_clusters = merged_lag[merged_lag['parent_cluster'] == merged_lag['child_cluster']]
    between_clusters = merged_lag[merged_lag['parent_cluster'] != merged_lag['child_cluster']]

    target_clusters = within_clusters

    grouped_by_cluster = target_clusters.groupby('parent_cluster')

    clusters = target_clusters['parent_cluster'].unique().tolist()
    cluster_summary = pd.DataFrame()

    pd_databases = ['community_agg_rank.pkl','RandomForest_agg_rank.pkl','Lasso_agg_rank.pkl']
    dict_databases = ['community_rank_data.pkl','RandomForest_rank_data.pkl','Lasso_rank_data.pkl']
    db = zip(pd_databases, dict_databases)
    agg_results = pd.DataFrame()
    #parsed_info = {}
    target_fp = 'marbach'
    """
    for pd_d, d_d in db:
        df = pd.read_pickle(pd_d)
        df = df[df['file_path'].str.contains(target_fp)]
        agg_results = agg_results.append(df)
        with open(d_d, 'rb') as infile:
            info = pickle.load(infile)
        
        parsed_info.update( {k: info[k] for k in df['result_path'].tolist() if k in info.keys()})
    """
    clusters.sort()
    current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
    output_file = "networks/community_marbach_clusters_" + current_time + ".tsv"
    output_file2 = "networks/SWING_community_marbach_clusters_" + current_time + ".tsv"
    output_file3 = "networks/SWING2_community_marbach_clusters_" + current_time + ".tsv"
    output_file4 = "networks/SWING3_community_marbach_clusters_" + current_time + ".tsv"
    #lagged_TFs = ['crp', 'fnr', 'fis', 'gadx', 'oxyr', 'iscr', 'csgd', 'torr']
    #lagged_TFs = ['crp','fnr', 'fis']
    #lagged_TFs = ['torr']
    #lagged_TFs = ['torr', 'cusr']
    lagged_TFs = None
    #for TF in lagged_TFs:
    current_group = merged_lag[(merged_lag['parent_cluster'] == 1) | (merged_lag['child_cluster'] == 1)]
    #current_group = current_group.append(merged_lag[(merged_lag['parent_cluster'] == 43) | (merged_lag['child_cluster'] == 43)])
    #current_group = merged_lag[merged_lag['parent'].isin(lagged_TFs)]
    #current_group_targets = current_group['child'].unique().tolist()
    #current_group_targets.extend(lagged_TFs)
    #current_group = merged_lag[merged_lag['child'].isin(current_group_targets)]
    #current_group = current_group[current_group['child'].isin(current_group_targets)]
    lagged_TFs = None
    clusterid = 25
    pr1,roc1 = run_subswing(current_group, td_window = 5, min_lag = 0, max_lag = 0, window_type = window_type, clusterid = clusterid, output_fn = output_file, tfs = lagged_TFs, alpha = .30, n_trees=500)
    pr2,roc2 = run_subswing(current_group, td_window = 4, min_lag = 1, max_lag = 1, window_type = window_type, clusterid = clusterid, output_fn = output_file2, tfs = lagged_TFs, alpha = .30, n_trees=500, pcs=4)
    #pr3,roc3 = run_subswing(current_group, td_window = 4, min_lag = 0, max_lag = 1, window_type = window_type, clusterid = clusterid, output_fn = output_file3, tfs = lagged_TFs)
    #pr4,roc4 = run_subswing(current_group, td_window = 3, min_lag = 1, max_lag = 2, window_type = window_type, clusterid = clusterid, output_fn = output_file4, tfs = lagged_TFs)
    """
    result5_list = []
    #alphas = np.arange(10,5000, 500)
    #alphas = np.arange(0.1,.4, 0.1)
    alphas = [1,2,3,4]
    for alpha in alphas:
        pr5,roc5 = run_subswing(current_group, td_window = 4, min_lag = 1, max_lag = 1, window_type = window_type, clusterid = clusterid, output_fn = output_file4, tfs = lagged_TFs, alpha=.29, pcs = alpha, n_trees = 1000)
        result5_list.append((alpha,pr5,roc5))
        print((alpha, pr5-pr1, roc5-roc1))
    pdb.set_trace()

    result5_list.sort(key=lambda tup:tup[2], reverse=True)
    """
    print('Diff pr: (negative is good):', pr1-pr2)
    print('diff roc:', roc1-roc2)
      #print(pr1,roc1,pr2,roc2,pr3,roc3,pr4,roc4,pr5,roc5,pr6,roc6)

    """                                                         
    subnet_info = extract_subnetwork(clusterid, merged_lag, parsed_info, agg_results)
    print('AUPR_Swing: ', subnet_info['swing_pr'])
    print('AUROC_Swing: ', subnet_info['swing_roc'])
    print('AUPR diff: ', subnet_info['baseline_pr'] - subnet_info['swing_pr'])
    print('AUROC diff: ', subnet_info['baseline_roc'] - subnet_info['swing_roc'])
    
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
                        #'swing_aupr4':pr5,
                        #'swing_auroc4':roc5,
                        #'swing_aupr5':pr6,
                        #'swing_auroc5':roc6

                        }
    cluster_summary = cluster_summary.append(cluster_result, ignore_index = True)
    """
    current_time = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')    
    cluster_summary.to_csv('cluster_summaries/sc_cluster_summary_within_community_c'+str(CLUSTER)+'_' + current_time + '.csv', header=True, index=False, sep='\t')


if __name__ == '__main__':
    if len(sys.argv) >= 2:
        window_type = str(sys.argv[1])
        CLUSTER = int(sys.argv[2])

    else:
        window_type = 'RandomForest'

    n_trials = 2
    for x in range(n_trials):
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
