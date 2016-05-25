import pdb
import pandas as pd
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
from collections import defaultdict


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


eco_go = parse_go()

clusters = pd.read_csv('../data/invitro/regulon_cluster_assignments7.csv',sep=',')


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


