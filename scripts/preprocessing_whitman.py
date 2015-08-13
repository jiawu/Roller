import sys,os,math
import pandas as pd
from PriorKnowledgeNetwork import PriorKnowledgeNetwork

import pdb


def check(col_names, df):
    for col_name in col_names:
        if col_name not in df.columns.values:
            raise ValueError(col_name + " does not exist in df")
def move_first(df):
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    return(df)

FOLDER_PATH = "../data/invitro/"
raw_whitfield_df = pd.read_csv(FOLDER_PATH+"whitfield_timeseries_raw.tsv",sep="\t")
whitfield_df = raw_whitfield_df.copy()

whitfield_map = pd.read_csv(FOLDER_PATH+"whitfield_1134_list.tsv",sep="\t")

parsed_genes = whitfield_df[whitfield_df['UID'].isin(whitfield_map['CLONEID'])]


if len(parsed_genes) != 1134:
    raise ValueError
parsed_genes = parsed_genes.append(whitfield_df[whitfield_df['UID']=='IMAGE:207808'])

parsed_genes['SYMBOL']=parsed_genes['NAME'].str.split(' ', expand=True, n=2).get(1)

#remove blank rows
#only process the genes with prior knowledge. Many of these transcripts describe genes that are not well known with no prior information. We will remove such transcripts from the dataset because there is no information about them in our gold standard.

parsed_genes = parsed_genes[parsed_genes['SYMBOL'] != '']
parsed_genes = parsed_genes[parsed_genes['SYMBOL'].notnull()]
## remove the blank columns

reference_pkn = PriorKnowledgeNetwork()
## remove genes with no prior knowledge
parsed_genes['BIOGRID'] = parsed_genes['SYMBOL'].apply(reference_pkn.get_candidate_ids)
parsed_genes['BIOGRID_length'] = parsed_genes['BIOGRID'].apply(len)

valid_genes = parsed_genes[parsed_genes['BIOGRID_length'] !=0]

#remove extra/uninformative columns
deleted_cols = ['UID','NAME','GWEIGHT','Fourier','Scaled Fourier', 'Sin', 'Cos', 'Atan2', 'Max Cor', 'AUTOCORR', 'EXP3 autocorr', 'EXP4 autocorr','Unnamed: 24', 'Unnamed: 51', 'Unnamed: 100', 'Unnamed: 120', 'BIOGRID','BIOGRID_length']

check(deleted_cols, valid_genes)

valid_genes.drop(deleted_cols,axis=1, inplace=True)

#many of the column names got messed up. So annoying. The timepoint 0.1 != 0.1 min. It is the second instance of "0". The next number 0.2 is the third instance of "0". Pandas renamed the columns to make them unique.

#create a separate time column, fix the pandas renaming timepoints issue
timepoints = valid_genes.columns.values.tolist()
del timepoints[-1]
timepoints = [math.floor(float(x)) for x in timepoints]

#transpose dataframe
valid_genes = valid_genes.set_index('SYMBOL').T
valid_genes.columns.name='Time'

valid_genes['Time'] = timepoints

#move the Time column first
cols = valid_genes.columns.tolist()
cols = cols[-1:] + cols[:-1]
valid_genes = valid_genes[cols]          

#create a perturbation column to label the experiments.
#first perturbation = timepoints[0:12]
perturbations = []
perturbations.extend(["Thy-Thy1" for x in range(0,12)])
#second perturbation= timepoints[12:12+26]
perturbations.extend(["Thy-Thy2" for x in range(12,(12+26))])
perturbations.extend(["Thy-Thy3" for x in range((12+26),(12+26+48))])
perturbations.extend(["Thy-Noc" for x in range((12+26+48),(12+26+48+19))])
perturbations.extend(["Shake" for x in range((12+26+48+19),(12+26+48+19+9))])


#remove all genes that have less than 114 timepoints
valid_genes = valid_genes.T.drop_duplicates().T
valid_genes['Perturbations'] = perturbations
valid_genes = valid_genes[valid_genes['Perturbations']=='Thy-Thy3']
valid_genes = valid_genes.dropna(axis=1,how='any')

#remove all genes except 9 target genes in Sambo, Shojaie, and Lozano
target_gene_list_muk = ['Time','CCNA2','E2F1','CCNB1','CCNE1','PCNA', 'CDC25A','CDC20','STK15','CKS2','PLK','DHFR','BUB1B','CCNF','BRCA1','TYMS','NFAT', 'CDC25C','CDC25B']

target_gene_list_shojaie = ['Time','CCNE1','CCNB1','CCNA2','RFC4','PCNA','E2F1','CDKN3','CDC6','CDC2']
check(target_gene_list_muk,valid_genes)
check(target_gene_list_shojaie,valid_genes)

muk_test = valid_genes[target_gene_list_muk]
shojaie_test = valid_genes[target_gene_list_shojaie]

#drop duplicates
muk_test = muk_test.T.groupby(level=0).first().T
shojaie_test = shojaie_test.T.groupby(level=0).first().T

#average time 0
muk_test=muk_test.groupby('Time').mean()
shojaie_test = shojaie_test.groupby('Time').mean()

muk_test['Time'] = [x for x in range(0,47)]
shojaie_test['Time'] = [x for x in range(0,47)]
#only for full dataset, not test sets. interpolate experiment  using cubic splines (as done in Lozano et al 2009). 
muk_test=move_first(muk_test)
shojaie_test=move_first(shojaie_test)

#save two raw data files

muk_test.to_csv(FOLDER_PATH+"whitfield_muk_timeseries.tsv",index=False,sep='\t')
shojaie_test.to_csv(FOLDER_PATH+"whitfield_shojaie_timeseries.tsv",index=False,sep='\t')

