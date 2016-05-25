import pandas as pd
import pdb
import scipy.stats as stats


def zscore_data(df):
    p = df.values
    z = pd.DataFrame(stats.zscore(p,axis=0,ddof=1),index=df.index, columns=df.columns)
    z['Time'] = df['Time']
    return(z)

db_path = '../data/invitro/net3_expression_data.tsv'
my_df = pd.read_csv(db_path, sep='\t')
my_df = my_df[~my_df['Time'].isnull()]

gp = my_df.groupby(['#Experiment','Time'])

#exp_list = [25,26,47,50,55,98, 105]

#my_df = my_df[my_df['#Experiment'].isin(exp_list)]

final_df = pd.DataFrame()

## Append certain rows with the same pertubation etc, alternating between repeats

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 25) & (my_df['Repeat'] == 1) ].iloc[:5])
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 25) & (my_df['Repeat'] == 2) ].iloc[:5])
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 26) & (my_df['Repeat'] == 1) ].iloc[:5])
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 26) & (my_df['Repeat'] == 2) ].iloc[:5])

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 47) & (my_df['Perturbations'].isnull())].iloc[:5])
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 47) & (my_df['Perturbations']=='P13')].iloc[:5])

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 50) & (my_df['Perturbations'] =='P8') & (my_df['Repeat'] == 1) ].iloc[:5] )
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 50) & (my_df['Perturbations'] =='P8') & (my_df['Repeat'] == 2) ].iloc[:5] )
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 50) & (my_df['Perturbations'] =='P8') & (my_df['Repeat'] == 3) ].iloc[:5] )

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 55) & (my_df['Perturbations'] =='P24') & (my_df['Repeat'] == 1) ].iloc[:5] )

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 98) & (my_df['DeletedGenes'].isnull()) & (my_df['Repeat'] == 1 )].iloc[:5] )

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 98) & (my_df['DeletedGenes'].isnull()) & (my_df['Repeat'] == 2 )].iloc[:5] )

final_df = final_df.append(my_df[ (my_df['#Experiment'] == 105)].iloc[:30] )
gene_names = pd.read_csv('../data/invitro/net3_gene_ids.tsv', sep='\t')

node_list = ['G%d'% (x) for x in range(1, 4512)]
node_list2 = gene_names['Name'].str.lower().tolist()
unmapped_df = final_df[['Time']+node_list]
unmapped_df.columns = ['Time'] + node_list2

om_df = pd.read_csv('../data/invitro/omranian_parsed_timeseries.tsv', sep='\t')
om_df = om_df[om_df['Time'] != 90]
intersecting_genes = set(om_df.columns.tolist()).intersection(set(unmapped_df.columns.tolist()))
intersecting_genes = sorted(list(intersecting_genes))
intersecting_genes.insert(0, intersecting_genes.pop(intersecting_genes.index('Time')))
mapped_df = unmapped_df[intersecting_genes]
norm_df = zscore_data(mapped_df)
# Change the time index so that it matches up with omranian...
x = [10,20,30,40,50]
t = [b for a in range(18) for b in x]
norm_df['Time'] = t

om_df_parsed = zscore_data(om_df[intersecting_genes])


om_df_parsed = om_df_parsed.append(norm_df)

om_df_parsed.to_csv('../data/invitro/omranian_parsed_timeseries.tsv', index=False, sep='\t')
pdb.set_trace()


