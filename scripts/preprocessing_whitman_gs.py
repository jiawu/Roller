from PriorKnowledgeNetwork import PriorKnowledgeNetwork
import pandas as pd
import networkx as nx
import itertools
import pdb
#import biogrid prior knowledge network

def generate_candidate_mappings(full_path):
    whitfield_df = pd.read_csv(full_path,sep="\t")
    items = whitfield_df.columns.values.tolist()
    items.pop(0)
    total_possible_edges = len(items) * (len(items)-1)

    #find list of candidate entrez IDs for each item
    overall_candidates = pd.DataFrame()
    for item in items:
        candidate_df = pd.DataFrame()
        candidate_df = candidate_df.append(reference_pkn.get_candidate_ids(item))
        candidate_df['Experiment Symbol']=item
        overall_candidates = overall_candidates.append(candidate_df)
    return(overall_candidates,items)
    
def generate_pkn(parsed_candidates, items):
    #using entrez ids generate a list of interactions using exhaustive search 
    # if the shortest path between two candidate ids is <X, then include it as an edge and stop looking for edges
    # make a list of possible edges
    edge_list = []
    entrez_dict = {}

    #iterate through each symbol, through each possible symbol
    #construct a dict to hold the experimental symbol, entrez values.

    for item in items:
        entrez_values = parsed_candidates[parsed_candidates['Experiment Symbol']==item]['Entrez'].tolist()
        entrez_dict[item] = entrez_values

    # two layers of iteration: (1) iterating through all node pairs, (2) within each pair, iterating through all entrez ID combos
    node_combinations = list(itertools.permutations(items,2))

    for source_symbol,target_symbol in node_combinations:
        id_combinations = list(itertools.product(entrez_dict[source_symbol],entrez_dict[target_symbol]))
        for source_id, target_id in id_combinations:
            try:
                path_len = reference_pkn.get_shortest_path(source=source_id,target=target_id)
                if path_len < 3:
                    edge_list.append((source_symbol,target_symbol))
                    break
                else:
                  continue
            except nx.NetworkXNoPath:
                continue

    pkn = nx.DiGraph()
    pkn.add_edges_from(edge_list)
    return(pkn)

def write_edge_list(outpath_name, pkn):
    formatted = nx.generate_edgelist(pkn, data=False, delimiter='\t')
    formatted = [s+"\t1" for s in formatted]

    with open(outpath_name,'w') as fout:
        fout.write("\n".join(formatted))  

reference_pkn = PriorKnowledgeNetwork()

#import data, determine the possible number of edges

FOLDER_PATH = "../data/invitro/"
full_path_shojaie = FOLDER_PATH+"whitfield_shojaie_timeseries.tsv"
full_path_muk = FOLDER_PATH+"whitfield_muk_timeseries.tsv"

overall_candidates_shojaie,shojaie_items = generate_candidate_mappings(full_path_shojaie)
overall_candidates_muk, muk_items = generate_candidate_mappings(full_path_muk)

#save entrez id map to file. parse manually (delete entrez IDs that do not match up). load entrez id map from file (manual parsing)
overall_candidates_shojaie.to_csv(FOLDER_PATH+"whitfield_shojaie_gs_candidate_mappings.tsv",sep="\t", index=False)

overall_candidates_muk.to_csv(FOLDER_PATH+"whitfield_muk_gs_candidate_mappings.tsv",sep="\t", index=False)

parsed_candidates_muk = pd.read_csv(FOLDER_PATH+"whitfield_muk_gs_corrected_mappings.tsv",sep="\t")
parsed_candidates_shojaie = pd.read_csv(FOLDER_PATH+"whitfield_shojaie_gs_corrected_mappings.tsv",sep="\t")

muk_pkn = generate_pkn(overall_candidates_muk,muk_items)
shojaie_pkn = generate_pkn(overall_candidates_shojaie,shojaie_items)

write_edge_list(FOLDER_PATH+"whitfield_muk_goldstandard.tsv",muk_pkn)
write_edge_list(FOLDER_PATH+"whitfield_shojaie_goldstandard.tsv",shojaie_pkn)

