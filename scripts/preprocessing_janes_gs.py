from PriorKnowledgeNetwork import PriorKnowledgeNetwork
import pandas as pd
import networkx as nx
import itertools
#import biogrid prior knowledge network
reference_pkn = PriorKnowledgeNetwork()

#import data, determine the possible number of edges

FOLDER_PATH = "../data/invitro/"
janes_df = pd.read_csv(FOLDER_PATH+"janes_timeseries_processed.txt",sep="\t")
items = janes_df.columns.values.tolist()
items.pop(0)
total_possible_edges = len(items) * (len(items)-1)

#find list of candidate entrez IDs for each item
overall_candidates = pd.DataFrame()
for item in items:
    candidate_df = pd.DataFrame()
    candidate_df = candidate_df.append(reference_pkn.get_candidate_ids(item))
    if "IRS1" in item:
        candidate_df = candidate_df.append(reference_pkn.get_candidate_ids('IRS1'))
    if "Caspase-8" in item:
        candidate_df = candidate_df.append(reference_pkn.get_candidate_ids('Caspase 8'))
        candidate_df = candidate_df.append(reference_pkn.get_candidate_ids('CASP8'))
    candidate_df['Experiment Symbol']=item
    overall_candidates = overall_candidates.append(candidate_df)
    
#save entrez id map to file. parse manually (delete entrez IDs that do not match up). load entrez id map from file (manual parsing)
overall_candidates.to_csv(FOLDER_PATH+"janes_gs_candidate_mappings.txt",sep="\t", index=False)

parsed_candidates = pd.read_csv(FOLDER_PATH+"janes_gs_corrected_mappings.txt",sep="\t")

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

formatted = nx.generate_edgelist(pkn, data=False, delimiter='\t')
formatted = [s+"\t1" for s in formatted]

with open(FOLDER_PATH+'janes_goldstandard.txt','w') as fout:
    fout.write("\n".join(formatted))  


