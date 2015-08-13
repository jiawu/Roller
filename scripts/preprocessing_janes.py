import os, sys, warnings
import pandas as pd

import pdb

def check(col_names, df):
    for col_name in col_names:
        if col_name not in df.columns.values:
            raise ValueError(col_name + " does not exist in df")
            
FOLDER_PATH = "../data/invitro/"
raw_janes_df = pd.read_csv(FOLDER_PATH+"janes_timeseries_raw.txt",sep="\t")
janes_df = raw_janes_df.copy()
## remove redundant data. see janes 2005, gaudet 2005, and ciaccio 2015.
## remove:
    ## tAkt: total Akt
    ## ptAkt: ratio of phosphorylated to total akt
    ## Akt: akt activity screened by kinase assay
    ## pAkt(IB): akt activity screened by immunoblot
    ## tEGFR: total EGFR
    ## ptEGFR : ratio of phospho to total EGFR
    ## ProC3: uncleaved casp3
    ## ProC8: uncleaved casp8

## remove apoptosis/viability data.
## remove:
    ## CC3/CCK
    ## Annexin
    ## PI
    ## Sub G1

deleted_columns = ['tAkt','ptAkt','Akt','pAkt (IB)','tEGFR','ptEGFR','ProC3','ProC8','CC3/CCK','Annexin','PI','Sub G1']

# check if deleted_columns is in df
check(deleted_columns, janes_df)
janes_df.drop(deleted_columns,axis=1, inplace=True)

# delete rows with NA datapoints
janes_df = janes_df.dropna()

# check number of replicates for each timepoint
grouped = janes_df.groupby('Time (min)')
grouped.count()

# combine the perturbation columns into one column
# there should be 11 unique perturbations
perturbations = [ 'TNF (ng/ml)', 'EGF (ng/ml)','Ins (ng/ml)', 'C225 (ug/ml)',
                  'IL-1ra (ug/ml)']
check(perturbations,janes_df)

janes_df['Perturbation']=janes_df.apply(lambda x: '%s_%s_%s_%s_%s' % (x[perturbations[0]],x[perturbations[1]],x[perturbations[2]],x[perturbations[3]],x[perturbations[4]]),axis=1)

# remove perturbation columns
janes_df.drop(perturbations,axis=1, inplace=True)
grouped = janes_df.groupby(['Time (min)','Perturbation']).count()

#unequal number of datapoints. time 0 has more 2x datapoints.
# average the replicates by: changing the replicates 4,5,6 to 1,2,3 and aggregation
janes_df['Replicate'] = janes_df['Replicate'].replace(4,1)
janes_df['Replicate'] = janes_df['Replicate'].replace(5,2)
janes_df['Replicate'] = janes_df['Replicate'].replace(6,3)

janes_df = janes_df.groupby(['Time (min)','Perturbation','Replicate'], as_index=False).mean()
grouped = janes_df.groupby('Time (min)')
grouped.count()
# remove replicate columns
janes_df.drop(['Replicate','Perturbation'],axis=1,inplace=True)

## make the column names nicer
    ## pAkt (AA) -> Akt
    ## ClvC8 -> Caspase-8
    ## pFKHR -> FKHR
    ## pMEK -> MEK
    ## pIRS1-636 -> IRS1-Ser
    ## pIRS1-896 -> IRS1-Tyr 
    # (I can't combine the two because the two series are somewhat uncorrelated with one another)
    ## pEGFR -> EGFR

renaming_key = {  'pAkt (AA)':'Akt',
                  'ClvC8':'Caspase-8',
                  'pFKHR':'FKHR',
                  'pMEK':'MEK',
                  'pIRS1-636':'IRS1-Ser',
                  'pIRS1-896':'IRS1-Tyr',
                  'pEGFR':'EGFR',
                  'Time (min)':'Time'}

check(list(renaming_key.keys()),janes_df)
janes_df.rename(columns=renaming_key, inplace=True)

janes_df.to_csv(FOLDER_PATH + 'janes_timeseries_processed.txt',sep='\t',index=False)

