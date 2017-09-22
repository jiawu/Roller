import pandas as pd
import pdb
import scipy.stats as stats
from scipy.interpolate import CubicSpline

def zscore_data(df):
    p = df.values
    z = pd.DataFrame(stats.zscore(p,axis=0,ddof=1),index=df.index, columns=df.columns)
    z['Time'] = df['Time']
    return(z)


def interpolate(df, target_timepoints):
    # interpolates between genes for 10 minute time steps
    # divides interpolated time series into # of non-overlapping timepoints
    # extracts target timepoints
    
    df['Time'] = (df['Time'] - df['Time'].iloc[0])*30
    cs = CubicSpline(df['Time'], df.iloc[:,8:])
    target_df = cs(target_timepoints)

    # create a new dataframe with:
    #pd.concat(df.iloc[0:6]*len(target_timepoints))
    i_gene = pd.DataFrame(target_df)
    i_gene.columns = df.iloc[:,8:].columns
    i_gene['Time'] = target_timepoints
    
    # fill in misc columns
    for col_name in df.columns:
        if col_name not in i_gene.columns:
            i_gene[col_name] = df[col_name].values[0]

    # return to original order
    i_gene = i_gene[df.columns]

    return(i_gene)
    

db_path = '../../data/invitro/net3_expression_data.tsv'
my_df = pd.read_csv(db_path, sep='\t')
my_df = my_df[~my_df['Time'].isnull()]

gp = my_df.groupby(['#Experiment','Time'])

#exp_list = [25,26,47,50,55,98, 105]

#my_df = my_df[my_df['#Experiment'].isin(exp_list)]

final_df = pd.DataFrame()

## Append certain rows with the same pertubation etc, alternating between repeats

#final_df = final_df.append(my_df[ (my_df['#Experiment'] == 25) & (my_df['Repeat'] == 1) ].iloc[:5])
#final_df = final_df.append(my_df[ (my_df['#Experiment'] == 25) & (my_df['Repeat'] == 2) ].iloc[:5])
#final_df = final_df.append(my_df[ (my_df['#Experiment'] == 26) & (my_df['Repeat'] == 1) ].iloc[:5])
#final_df = final_df.append(my_df[ (my_df['#Experiment'] == 26) & (my_df['Repeat'] == 2) ].iloc[:5])
ts = my_df[ (my_df['#Experiment'] == 47) & (my_df['Perturbations'].isnull())]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))
final_df = final_df.append(interpolate(ts, [50,60,70,80,90]))
final_df = final_df.append(interpolate(ts, [100,110,120,130,140]))
final_df = final_df.append(interpolate(ts, [150,160,170,180,190]))

ts = my_df[ (my_df['#Experiment'] == 47) & (my_df['Perturbations']=='P13')]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))
final_df = final_df.append(interpolate(ts, [50,60,70,80,90]))
final_df = final_df.append(interpolate(ts, [100,110,120,130,140]))
final_df = final_df.append(interpolate(ts, [150,160,170,180,190]))

ts = my_df[ (my_df['#Experiment'] == 49) & (my_df['Perturbations'].isnull())]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))
final_df = final_df.append(interpolate(ts, [50,60,70,80,90]))


# first time point is shared between the other two experiments
temp_t0 = my_df[ (my_df['#Experiment'] == 49) & (my_df['Perturbations'].isnull())].iloc[0,:]

ts = pd.DataFrame().append(temp_t0)
ts = ts.append(my_df[ (my_df['#Experiment'] == 49) & (my_df['Perturbations'] == 'P16')].iloc[:4])
ts = ts[final_df.columns]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))
final_df = final_df.append(interpolate(ts, [50,60,70,80,90]))

ts = pd.DataFrame().append(temp_t0)
ts = ts.append(my_df[ (my_df['#Experiment'] == 49) & (my_df['Perturbations'] == 'P17')].iloc[:4])
ts = ts[final_df.columns]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))
final_df = final_df.append(interpolate(ts, [50,60,70,80,90]))


final_df = final_df.append(my_df[ (my_df['#Experiment'] == 50) & (my_df['Perturbations'] =='P8') & (my_df['Repeat'] == 1) ].iloc[:5] )
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 50) & (my_df['Perturbations'] =='P8') & (my_df['Repeat'] == 2) ].iloc[:5] )
final_df = final_df.append(my_df[ (my_df['#Experiment'] == 50) & (my_df['Perturbations'] =='P8') & (my_df['Repeat'] == 3) ].iloc[:5] )

#final_df = final_df.append(my_df[ (my_df['#Experiment'] == 55) & (my_df['Perturbations'] =='P24') & (my_df['Repeat'] == 1) ].iloc[:5] )

#final_df = final_df.append(my_df[ (my_df['#Experiment'] == 98) & (my_df['DeletedGenes'].isnull()) & (my_df['Repeat'] == 1 )].iloc[:5] )

#final_df = final_df.append(my_df[ (my_df['#Experiment'] == 98) & (my_df['DeletedGenes'].isnull()) & (my_df['Repeat'] == 2 )].iloc[:5] )

ts= my_df[ (my_df['#Experiment'] == 105)].iloc[0:5]
ts = ts[final_df.columns]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))

ts= my_df[ (my_df['#Experiment'] == 105)].iloc[5:10]
ts = ts[final_df.columns]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))

ts= my_df[ (my_df['#Experiment'] == 105)].iloc[10:15]
ts = ts[final_df.columns]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))

ts= my_df[ (my_df['#Experiment'] == 105)].iloc[15:20]
ts = ts[final_df.columns]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))

ts= my_df[ (my_df['#Experiment'] == 105)].iloc[20:25]
ts = ts[final_df.columns]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))

ts= my_df[ (my_df['#Experiment'] == 105)].iloc[25:30]
ts = ts[final_df.columns]
final_df = final_df.append(interpolate(ts, [0,10,20,30,40]))

gene_names = pd.read_csv('../../data/invitro/net3_gene_ids.tsv', sep='\t')

node_list = ['G%d'% (x) for x in range(1, 4512)]
node_list2 = gene_names['Name'].str.lower().tolist()
unmapped_df = final_df[['Time']+node_list]
unmapped_df.columns = ['Time'] + node_list2

om_df = pd.read_csv('../../data/invitro/iomranian_parsed_timeseries.tsv', sep='\t')
om_df = om_df[om_df['Time'] != 90]
intersecting_genes = set(om_df.columns.tolist()).intersection(set(unmapped_df.columns.tolist()))
intersecting_genes = sorted(list(intersecting_genes))
intersecting_genes.insert(0, intersecting_genes.pop(intersecting_genes.index('Time')))


mapped_df = unmapped_df[intersecting_genes]
norm_df = zscore_data(mapped_df)

# Change the time index so that it matches up with omranian...
x = [10,20,30,40,50]
t = [b for a in range(int(norm_df.shape[0]/5)) for b in x]
norm_df['Time'] = t

om_df_parsed = zscore_data(om_df[intersecting_genes])


om_df_parsed = om_df_parsed.append(norm_df)

om_df_parsed.to_csv('../../data/invitro/iomranian_parsed_timeseries.tsv', index=False, sep='\t')


