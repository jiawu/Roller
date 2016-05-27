import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import interp1d

import sys
from datetime import datetime
import numpy as np

sys.path.append("../pipelines")
import Pipelines as tdw
import Swing.util.lag_identification as lag_id
from Swing.util.Evaluator import Evaluator

import pdb

def create_df(raw_data_path):
    raw_data = pd.read_csv(raw_data_path, sep = '\t')
    time = [10,20,30,40,50]

    raw_data = raw_data.transpose()
    probe_names=raw_data.iloc[0]
    raw_data = raw_data[1:][:]

    biol_reps = 3

    final_frame = pd.DataFrame()

    for i in range(3):
        rep_frame = raw_data[i::3][:].iloc[:5]
        final_frame = final_frame.append(rep_frame)

#make time series
    final_time_series = time + time + time
    final_frame['Time'] = final_time_series

    cols = final_frame.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    final_frame=final_frame[cols]

#mapping the names to the map and the gold standard
    name_map = pd.read_csv('../data/invitro/ecoli_gene_map.txt', sep = '\t')

#map the probe name to the gene name
    gene_list = []
    for name in probe_names:
        gene_name = name_map[name_map['ID'] == name]['ORF']
        gene_list.append(gene_name.values[0].lower())

    gene_list.insert(0, 'Time')
    final_frame.columns = gene_list

    return(final_frame)

#check if gene name exists in map

#raw_data_list = ['../data/invitro/omranian_coldstress.txt','../data/invitro/omranian_heatstress.txt','../data/invitro/omranian_control.txt','../data/invitro/omranian_oxidativestress.txt']
raw_data_list = ['../data/invitro/omranian_control.txt','../data/invitro/omranian_coldstress.txt','../data/invitro/omranian_heatstress.txt', '../data/invitro/omranian_oxidativestress.txt']
#raw_data_list = ['../data/invitro/omranian_oxidativestress.txt']

overall_df = pd.DataFrame()
"""
for raw_data in raw_data_list:
    df = create_df(raw_data)
    overall_df = overall_df.append(df)
"""
overall_df = pd.read_csv('../data/invitro/omranian_parsed_timeseries.tsv',sep='\t')

genes = overall_df.columns[1:].tolist()

with open('../data/invitro/omranian_tf_list.tsv','r') as f:
    tf_list = f.read().splitlines()
with open('../data/invitro/omranian_target_list.tsv','r') as f:
    target_list = f.read().splitlines()

gs_list = tf_list+target_list

gs_list = list(set(gs_list))


not_measured_genes = list(set(gs_list)-set(genes))
#remove these genes from gold standard

parsed_gs = pd.read_csv('../data/invitro/omranian_parsed_goldstandard.tsv',sep='\t', header=None)

parsed_gs.columns = ['regulator', 'target', 'effect']
parsed_gs = parsed_gs[~parsed_gs['regulator'].isin(not_measured_genes)]
parsed_gs = parsed_gs[~parsed_gs['target'].isin(not_measured_genes)]

parsed_gs = parsed_gs[~(parsed_gs['regulator']==parsed_gs['target'])]

parsed_gs.to_csv('../data/invitro/omranian_parsed_goldstandard.tsv',sep='\t',index=False, header=False)

parsed_gs_list = list(set(parsed_gs['regulator'].tolist()+parsed_gs['target'].tolist()))

with open('../data/invitro/omranian_all_genes_list.tsv', 'w') as outfile:
    outfile.write("\n".join(parsed_gs_list))

with open('../data/invitro/omranian_parsed_tf_list.tsv', 'w') as outfile:
    outfile.write("\n".join(parsed_gs['regulator'].unique().tolist()))

not_in_gs = list(set(genes) - set(gs_list))
#remove these genes from overall_df

in_df_in_gs = ['Time'] + list(set(gs_list).intersection(genes))
final_overall_df = overall_df[in_df_in_gs]

final_overall_df.to_csv('../data/invitro/omranian_parsed_timeseries.tsv', index=False,sep='\t')

my_eval = Evaluator('../data/invitro/omranian_parsed_goldstandard.tsv')
