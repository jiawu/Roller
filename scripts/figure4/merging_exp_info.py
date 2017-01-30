import pdb
import pandas as pd


db_path = '../data/invitro/net4_chip_features.tsv'
my_df = pd.read_csv(db_path, sep='\t')
db_path2 = '../data/invitro/net4_expression_data.tsv'
my_df2 = pd.read_csv(db_path2, sep='\t')
my_df = my_df.join(my_df2)

master_map = pd.read_csv('../data/invitro/marbach_gene_ids.tsv', sep='\t')
master_map.columns = ['anonID', 'geneid']
map_dict = master_map.set_index('anonID').T.to_dict('record')

# replace gs file with IDs
parsed_df = my_df.rename(columns=map_dict[0])
main_columns = parsed_df.columns[8:]
main_columns = [x for x in main_columns if 'decoy' not in x]

key_df = pd.read_csv('../data/invitro/yeast_exp_key.tsv',sep='\t')
key_df=key_df.rename(columns = {'Unnamed: 0':'author', 'Unnamed: 1': 'exp_info'})
key_df['hash'] = key_df[main_columns].astype(str).sum(axis=1)
parsed_df['hash'] = parsed_df[main_columns].astype(str).sum(axis=1)

exp_map = pd.merge(parsed_df, key_df, on=['hash'], how='inner')
exp_map = exp_map[['#Experiment','author','exp_info']]
exp_map.to_csv('Yeast_exp_map.tsv', sep='\t',index=False)
