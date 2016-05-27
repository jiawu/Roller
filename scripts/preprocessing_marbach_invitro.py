import pandas as pd
import pdb
import scipy.stats as stats


def zscore_data(df):
    p = df.values
    z = pd.DataFrame(stats.zscore(p,axis=0,ddof=1),index=df.index, columns=df.columns)
    z['Time'] = df['Time']
    return(z)

db_path = '../data/invitro/net4_chip_features.tsv'
my_df = pd.read_csv(db_path, sep='\t')
db_path2 = '../data/invitro/net4_expression_data.tsv'
my_df2 = pd.read_csv(db_path2, sep='\t')
my_df = my_df.join(my_df2)
my_df = my_df[~my_df['Time'].isnull()]

gp = my_df.groupby(['#Experiment','Time'])

#exp_list = [25,26,47,50,55,98, 105]

#my_df = my_df[my_df['#Experiment'].isin(exp_list)]

final_df = pd.DataFrame()

## Append certain rows with the same pertubation etc, alternating between repeats

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 26) & (my_df['Repeat'] == 1) ].iloc[:6])
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 30) & (my_df['Repeat'] == 1) ].iloc[:6])
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 48) & (my_df['Repeat'] == 1) ].iloc[:12])
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 49) & (my_df['Repeat'] == 1) &(my_df['Perturbations'].isnull())].iloc[:6])
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 49) & (my_df['Repeat'] == 1) &(my_df['Perturbations']=='P19')].iloc[:36])

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 51) & (my_df['Repeat'] == 1) &(my_df['Perturbations'].isnull())].iloc[:6])
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 51) & (my_df['Repeat'] == 1) &(my_df['Perturbations']=='P24')].iloc[:6])

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 58) & (my_df['Repeat'] == 1)&(my_df['DeletedGenes'].isnull())].iloc[:30])

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 58) & (my_df['Repeat'] == 1)&(my_df['DeletedGenes']=='G3606')].iloc[:30])


master_map = pd.read_csv('../data/invitro/marbach_gene_ids.tsv', sep='\t')
master_map.columns = ['anonID', 'geneid']
map_dict = master_map.set_index('anonID').T.to_dict('record')
# replace gs file with IDs
parsed_df = final_df.rename(columns=map_dict[0])
# get time column and all the genes
col_names = ['Time'] + master_map['geneid'].values.tolist()
parsed_df = parsed_df[col_names]

gs = pd.read_csv('../data/invitro/marbach_parsed_goldstandard.tsv',sep='\t', header=None)

all_genes = gs.iloc[:,0].unique().tolist() + gs.iloc[:,1].unique().tolist()
all_genes_gs = set(all_genes)
all_genes_ts = set(master_map['geneid'].values.tolist())

#get intersection
shared_genes = all_genes_ts.intersection(all_genes_gs)
col_names = ['Time'] + list(shared_genes)
parsed_df = parsed_df[col_names]


tf_list = pd.read_csv('../data/invitro/marbach_tf_list.tsv',sep='\t', header=None)
shared_tfs = list(shared_genes.intersection(set(tf_list.iloc[:,0].tolist())))

with open('../data/invitro/marbach_parsed_tf_list.tsv', 'w') as outfile:
    outfile.write("\n".join(shared_tfs))
with open('../data/invitro/marbach_all_genes_list.tsv', 'w') as outfile:
    outfile.write("\n".join(list(shared_genes)))


# zscore the data
# check if gold standard has all the names, or gold standard is measuring these species... remove the decoys for example
# parse the tf list to have the proper mappings, names

# marbach_parsed_goldstandard.tsv x
# marbach_parsed_timeseries.tsv 
# marbach_parsed_tf_list.tsv x
# marbach_all_gene_list.tsv x

# marbach_signed_parsed_goldstandard.tsv

norm_df = zscore_data(parsed_df)
# Change the time index so that it matches up with omranian...
x = [10,20,30,40,50,60]
t = [b for a in range(23) for b in x]
norm_df['Time'] = t


norm_df.to_csv('../data/invitro/marbach_parsed_timeseries.tsv', index=False, sep='\t')


