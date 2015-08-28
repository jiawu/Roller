import os, sys, warnings
import pandas as pd
import networkx as nx
import numpy as np
import pdb

"""Generates a prior knowledge network in the form of a networkx model"""

class PriorKnowledgeNetwork:

    def __init__(self):
        FOLDER_PATH = "../data/invitro/"
        raw_biogrid_df = pd.read_csv(FOLDER_PATH+"biogrid_full.tsv",sep="\t")
        self.bg_df = raw_biogrid_df.copy()
        self.bg_entrez = raw_biogrid_df.copy()
        self.bg_entrez = self.get_entrez_mapping_df(self.bg_entrez)

        #note that these are often bait/prey interactions and not functional interactions
        bait_list = self.bg_df['Entrez Gene Interactor A']
        prey_list = self.bg_df['Entrez Gene Interactor B']

        interaction_list = list(zip(bait_list,prey_list))

        self.directed_graph = nx.DiGraph()
        self.directed_graph.add_edges_from(interaction_list)

    def get_entrez_mapping_df(self,bg_df):
        # get the following columns: entrez gene interactor a, b, official symbol a,b, synonyms a,b
        # create a dataframe that looks like (entrez, official symbol, synonyms)
        # remove duplicates
        A_interactions= bg_df[['Entrez Gene Interactor A','Official Symbol Interactor A', 'Synonyms Interactor A']]
        B_interactions= bg_df[['Entrez Gene Interactor B','Official Symbol Interactor B', 'Synonyms Interactor B']]
        header = ['Entrez','Official Symbol','Synonyms'] 
        A_interactions.columns = header
        B_interactions.columns = header

        all_interactions = A_interactions.append(B_interactions)
        all_interactions.drop_duplicates(inplace=True)
        return(all_interactions)

    def get_shortest_path(self,source=None, target=None):
        path = nx.shortest_path(self.directed_graph,source=source, target=target)
        return(len(path))
    
    def get_candidate_ids(self, gene_name):
        search_col = ['Official Symbol', 'Synonyms']
        mydf = self.bg_entrez
        mask = np.column_stack([mydf[col].str.contains("(?i)"+gene_name, na=False) for col in mydf if col in search_col])
        candidate_gene_table = mydf.loc[mask.any(axis=1)]
        return(candidate_gene_table)

    def check(self, col_names, df):
        for col_name in col_names:
            if col_name not in df.columns.values:
                raise ValueError(col_name + " does not exist in df")
                

