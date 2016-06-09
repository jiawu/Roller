import pdb
import pandas as pd
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
from collections import defaultdict


def parse_go():
    go = pd.read_csv('../data/invitro/gene_association.sgd', sep='\t', comment='!')
    master_map = pd.read_csv('../data/invitro/marbach_gene_ids.tsv', sep='\t')
    gs = pd.read_csv('../data/invitro/marbach_goldstandard.tsv',sep='\t')
    tfs = pd.read_csv('../data/invitro/marbach_tf_list.tsv',sep='\t')


    gs.columns = ['regulator', 'target', 'exists']

    genes_go = go.iloc[:,[2,4,10]]
    genes_go.columns = ['name','GO_ID', 'alias']
    #parse alias
    genes_go['parsed_alias']= genes_go['alias'].str.split('|').str[0]
    
    # first make the eco go annotations
    genes = genes_go['parsed_alias'].tolist()
    go_id = genes_go['GO_ID'].tolist()
    go_tuple = list(zip(genes,go_id))
    eco_go = defaultdict(list)
    for genes,go_id in go_tuple:
        eco_go[genes].append(go_id)

    # next, make the master table that maps geneid to anonymized gene name to real genename in eco_go
    master_map.columns = ['anonID', 'geneid']
    map_dict = master_map.set_index('anonID').T.to_dict()
    # replace gs file with IDs
    parsed_gs = gs.replace(map_dict)
    parsed_gs.to_csv('../data/invitro/marbach_parsed_goldstandard.tsv',sep='\t',header=None,index=False)

    return(eco_go)

def get_clist(clusterid, cluster_table):
    clist = cluster_table[cluster_table['__glayCluster'] == clusterid]['name'].tolist()
    return(clist)


eco_go = parse_go()
pdb.set_trace()

clusters = pd.read_csv('../data/invitro/yeast_cluster_assign5.csv',sep=',')


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
r = []
for clusterid in valid_clusters:
    genes_0 = get_clist(clusterid, clusters)
    goea_results_all = goeaobj.run_study(genes_0)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_by < 0.05]
    print(len(goea_results_sig))
    ratios = []
    for result in goea_results_sig:
        result_tup = (clusterid,result.name, result.ratio_in_study)
        ratios.append(result_tup)

    ratios = sorted(ratios, key=lambda tup: tup[1][2][0])
    my_result = { 'clusterid': clusterid,
                  'ratios': ratios,
                  'n_genes': len(genes_0)}

    r.append(my_result)

count_ontology = 0
count_no = 0
for result in r:
    # get number of modules over size 10, and has an ontology
    # get number of modules over size 10 with no ontology
    if result['n_genes'] > 9:
        if not result['ratios']:
            count_no += 1
        else:
            count_ontology +=1
            print(result)

print("# with ont", count_ontology)
print("# w/o ont", count_no)


